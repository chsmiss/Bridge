#!/usr/bin/env python3
"""Score GFA junctions with a small, transparent nucleotide language model.

This is a reproducible adapter baseline, not a replacement for GENERator/CARBon.  It
trains an order-k Markov model on graph segments and measures whether the first bases of
a target are more probable when conditioned on the source tail than when the target is
started independently.  Neural DNA models should emit the same source/target/
dna_lm_delta TSV schema.
"""

from __future__ import annotations

import argparse
import collections
import math
import sys
from pathlib import Path
from typing import DefaultDict, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "ACGT"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def parse_overlap(value: str) -> int:
    if value == "*":
        return 0
    if value.endswith("M") and value[:-1].isdigit():
        return int(value[:-1])
    return 0


def read_gfa(path: Path) -> Tuple[Dict[str, str], List[Tuple[str, str, int]]]:
    segments: Dict[str, str] = {}
    links: List[Tuple[str, str, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("H\t"):
                continue
            fields = line.split("\t")
            if fields[0] == "S":
                if len(fields) < 3 or fields[2] == "*":
                    raise ValueError(f"{path}:{line_number}: sequence-bearing S line required")
                segments[fields[1]] = fields[2].upper()
            elif fields[0] == "L" and len(fields) >= 6:
                if fields[2] == "+" and fields[4] == "+":
                    links.append((fields[1], fields[3], parse_overlap(fields[5])))
    return segments, links


class MarkovModel:
    def __init__(self, order: int, alpha: float) -> None:
        self.order = order
        self.alpha = alpha
        self.counts: DefaultDict[str, DefaultDict[str, int]] = collections.defaultdict(
            lambda: collections.defaultdict(int)
        )
        self.totals: DefaultDict[str, int] = collections.defaultdict(int)

    def observe(self, sequence: str) -> None:
        sequence = sequence.upper()
        padded = "^" * self.order + sequence
        for index, base in enumerate(sequence):
            if base not in BASES:
                continue
            context = padded[index : index + self.order]
            if any(symbol not in BASES + "^" for symbol in context):
                continue
            self.counts[context][base] += 1
            self.totals[context] += 1

    def log_probability(self, base: str, context: str) -> float:
        if base not in BASES:
            return math.log(0.25)
        context = context[-self.order :].rjust(self.order, "^")
        total = self.totals.get(context, 0)
        count = self.counts.get(context, {}).get(base, 0)
        return math.log((count + self.alpha) / (total + self.alpha * 4.0))

    def score_extension(self, prefix: str, extension: str, limit: int) -> float:
        context = prefix[-self.order :].rjust(self.order, "^")
        score = 0.0
        used = 0
        for base in extension[:limit].upper():
            if base not in BASES:
                context = (context + "N")[-self.order :]
                continue
            score += self.log_probability(base, context)
            context = (context + base)[-self.order :]
            used += 1
        return score / max(1, used)


def train_model(segments: Mapping[str, str], order: int, alpha: float) -> MarkovModel:
    model = MarkovModel(order=order, alpha=alpha)
    for sequence in segments.values():
        model.observe(sequence)
        model.observe(reverse_complement(sequence))
    return model


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gfa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--order", type=int, default=5)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--score-nt", type=int, default=96)
    parser.add_argument("--clip", type=float, default=3.0, help="clip per-base delta before output")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if not 1 <= args.order <= 12:
            raise ValueError("--order must be in [1, 12]")
        if args.alpha <= 0.0:
            raise ValueError("--alpha must be positive")
        if args.score_nt < 1:
            raise ValueError("--score-nt must be positive")
        segments, links = read_gfa(args.gfa)
        if not segments:
            raise ValueError("GFA contains no segments")
        model = train_model(segments, order=args.order, alpha=args.alpha)
        with args.output.open("w", encoding="utf-8") as handle:
            handle.write(
                "source\ttarget\tdna_lm_delta\tjoined_logp_per_base\tindependent_logp_per_base\tmodel\n"
            )
            for source, target, overlap in links:
                if source not in segments or target not in segments:
                    continue
                left = segments[source]
                right = segments[target][min(overlap, len(segments[target])) :]
                joined = model.score_extension(left, right, args.score_nt)
                independent = model.score_extension("", right, args.score_nt)
                delta = max(-args.clip, min(args.clip, joined - independent))
                handle.write(
                    f"{source}\t{target}\t{delta:.6f}\t{joined:.6f}\t"
                    f"{independent:.6f}\tmarkov_order_{args.order}\n"
                )
        print(
            f"DNA junction baseline: segments={len(segments)} links={len(links)} order={args.order}",
            file=sys.stderr,
        )
        return 0
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
