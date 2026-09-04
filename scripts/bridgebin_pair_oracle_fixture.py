#!/usr/bin/env python3
"""Generate a sparse truth-derived pair-score fixture for testing BridgeBin v2.1 plumbing.

THIS IS NOT A BIOLOGICAL MODEL AND MUST NOT BE USED AS A PERFORMANCE RESULT.
It exists only to answer a software question: if a future DNA/GENERanno/ESM-C pair head
can confidently distinguish two coverage-identical genomes, does the Rust refinement
layer actually split the mixed bin without destroying pure bins?

For each truth genome, long contigs are selected as anchors. Every eligible contig gets
high-p_same links to same-genome anchors and low-p_same links to a small number of anchors
from every other genome. The output schema is exactly the production ``--pair-scores``
interface.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--anchors-per-genome", type=int, default=4)
    parser.add_argument("--negative-anchors-per-genome", type=int, default=2)
    parser.add_argument("--positive", type=float, default=0.995)
    parser.add_argument("--negative", type=float, default=0.005)
    parser.add_argument("--confidence", type=float, default=0.999)
    return parser.parse_args(argv)


def truth_rows(path: Path) -> Iterable[Tuple[str, str, int, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            genome = (row.get("genome") or row.get("species") or "").strip()
            if not contig or not genome:
                continue
            length = int(row.get("length") or 0)
            eligible = int(row.get("eligible") or 1)
            yield contig, genome, length, eligible


def ordered(left: str, right: str) -> Tuple[str, str]:
    return (left, right) if left <= right else (right, left)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.anchors_per_genome < 1 or args.negative_anchors_per_genome < 1:
        raise SystemExit("anchor counts must be positive")
    if not (0.0 <= args.negative < args.positive <= 1.0):
        raise SystemExit("require 0 <= negative < positive <= 1")
    if not 0.0 <= args.confidence <= 1.0:
        raise SystemExit("confidence must be in [0,1]")

    by_genome: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    for contig, genome, length, eligible in truth_rows(args.truth):
        if eligible:
            by_genome[genome].append((contig, length))
    if len(by_genome) < 2:
        raise SystemExit("fixture requires at least two truth genomes")
    for values in by_genome.values():
        values.sort(key=lambda item: (-item[1], item[0]))

    anchors = {
        genome: [contig for contig, _ in values[: args.anchors_per_genome]]
        for genome, values in by_genome.items()
    }
    negatives = {
        genome: values[: args.negative_anchors_per_genome]
        for genome, values in anchors.items()
    }

    rows: Dict[Tuple[str, str], Tuple[float, str]] = {}
    for genome, values in sorted(by_genome.items()):
        same_anchors = anchors[genome]
        for contig, _length in values:
            for anchor in same_anchors:
                if contig == anchor:
                    continue
                rows[ordered(contig, anchor)] = (args.positive, "oracle_same_anchor")
            for other_genome, other_anchors in sorted(negatives.items()):
                if other_genome == genome:
                    continue
                for anchor in other_anchors:
                    if contig == anchor:
                        continue
                    rows[ordered(contig, anchor)] = (args.negative, "oracle_cross_anchor")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["left", "right", "p_same", "confidence", "model"])
        for (left, right), (probability, source) in sorted(rows.items()):
            writer.writerow(
                [
                    left,
                    right,
                    f"{probability:.6f}",
                    f"{args.confidence:.6f}",
                    f"oracle_fixture:{source}",
                ]
            )
    print(
        f"bridgebin-oracle-pairs: genomes={len(by_genome)} pairs={len(rows)} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
