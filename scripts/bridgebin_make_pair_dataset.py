#!/usr/bin/env python3
"""Build labelled same-genome training pairs with deliberate hard-negative mining.

The dataset is sparse and anchor-based: each eligible contig is paired with several
same-genome anchors and with the most coverage-compatible anchors from other genomes.
This deliberately over-samples the exact failure mode that hurts metagenomic binning:
different genomes with nearly indistinguishable abundance profiles.

Truth is used only for labels/group metadata and for choosing positive/negative genomes;
it is never used as a model feature. Output includes left_genome/right_genome so
``bridgebin_pair_head.py`` can hold out complete genomes during validation.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Record = Tuple[str, str, int]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coverage", type=Path)
    parser.add_argument("--anchors-per-genome", type=int, default=24)
    parser.add_argument("--positive-per-contig", type=int, default=4)
    parser.add_argument("--negative-genomes-per-contig", type=int, default=4)
    parser.add_argument("--negative-anchors-per-genome", type=int, default=2)
    parser.add_argument("--min-length", type=int, default=1500)
    parser.add_argument("--max-contigs-per-genome", type=int, default=0)
    parser.add_argument("--seed", type=int, default=43)
    return parser.parse_args(argv)


def read_truth(path: Path, min_length: int) -> List[Record]:
    out: List[Record] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            genome = (row.get("genome") or row.get("species") or row.get("bin") or "").strip()
            if not contig or not genome:
                continue
            length = int(float(row.get("length") or 0))
            eligible_raw = (row.get("eligible") or "1").strip().lower()
            eligible = eligible_raw not in {"0", "false", "no"}
            if eligible and length >= min_length:
                out.append((contig, genome, length))
    return out


def read_coverage(path: Optional[Path]) -> Dict[str, List[float]]:
    if path is None:
        return {}
    result: Dict[str, List[float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError(f"{path}: coverage table needs contig and sample columns")
        for fields in reader:
            if len(fields) < 2 or not fields[0].strip():
                continue
            values = [float(value) for value in fields[1:]]
            if all(math.isfinite(value) and value >= 0.0 for value in values):
                result[fields[0].strip()] = values
    return result


def coverage_similarity(left: Sequence[float], right: Sequence[float]) -> Optional[float]:
    if not left or len(left) != len(right):
        return None
    distance = sum(abs(math.log((a + 0.5) / (b + 0.5))) for a, b in zip(left, right)) / len(left)
    return math.exp(-distance / 0.80)


def ordered(left: str, right: str) -> Tuple[str, str]:
    return (left, right) if left <= right else (right, left)


def cap_genomes(records: List[Record], max_per_genome: int, seed: int) -> List[Record]:
    if max_per_genome <= 0:
        return records
    grouped: Dict[str, List[Record]] = defaultdict(list)
    for record in records:
        grouped[record[1]].append(record)
    kept: List[Record] = []
    for genome, values in sorted(grouped.items()):
        rng = random.Random(f"{seed}:{genome}")
        values = sorted(values, key=lambda record: (-record[2], record[0]))
        if len(values) > max_per_genome:
            # Preserve long anchors but sample the remaining diversity deterministically.
            head = values[: min(max_per_genome // 3, len(values))]
            tail = values[len(head) :]
            rng.shuffle(tail)
            values = head + tail[: max_per_genome - len(head)]
        kept.extend(values)
    return kept


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    for name in (
        "anchors_per_genome",
        "positive_per_contig",
        "negative_genomes_per_contig",
        "negative_anchors_per_genome",
    ):
        if getattr(args, name) < 1:
            raise SystemExit(f"--{name.replace('_', '-')} must be positive")

    records = cap_genomes(read_truth(args.truth, args.min_length), args.max_contigs_per_genome, args.seed)
    coverage = read_coverage(args.coverage)
    by_genome: Dict[str, List[Record]] = defaultdict(list)
    genome_of: Dict[str, str] = {}
    for record in records:
        by_genome[record[1]].append(record)
        genome_of[record[0]] = record[1]
    if len(by_genome) < 3:
        raise SystemExit("need at least three genomes for leakage-aware pair training")
    for values in by_genome.values():
        values.sort(key=lambda record: (-record[2], record[0]))

    anchors: Dict[str, List[str]] = {
        genome: [record[0] for record in values[: args.anchors_per_genome]]
        for genome, values in by_genome.items()
    }
    pairs: Dict[Tuple[str, str], Tuple[int, str, str, Optional[float], str]] = {}
    rng = random.Random(args.seed)

    for contig, genome, _length in sorted(records):
        same = [anchor for anchor in anchors[genome] if anchor != contig]
        rng.shuffle(same)
        for anchor in same[: args.positive_per_contig]:
            key = ordered(contig, anchor)
            similarity = coverage_similarity(coverage.get(contig, []), coverage.get(anchor, []))
            pairs[key] = (1, genome_of[key[0]], genome_of[key[1]], similarity, "positive_anchor")

        genome_scores: List[Tuple[float, str, List[Tuple[float, str]]]] = []
        for other_genome, other_anchors in anchors.items():
            if other_genome == genome:
                continue
            ranked: List[Tuple[float, str]] = []
            for anchor in other_anchors:
                similarity = coverage_similarity(coverage.get(contig, []), coverage.get(anchor, []))
                ranked.append((similarity if similarity is not None else 0.0, anchor))
            ranked.sort(key=lambda item: (-item[0], item[1]))
            best = ranked[0][0] if ranked else 0.0
            genome_scores.append((best, other_genome, ranked))
        genome_scores.sort(key=lambda item: (-item[0], item[1]))

        for _best, _other_genome, ranked in genome_scores[: args.negative_genomes_per_contig]:
            for similarity, anchor in ranked[: args.negative_anchors_per_genome]:
                key = ordered(contig, anchor)
                pairs[key] = (
                    0,
                    genome_of[key[0]],
                    genome_of[key[1]],
                    similarity if coverage else None,
                    "coverage_hard_negative" if coverage else "negative_anchor",
                )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "left",
                "right",
                "label",
                "left_genome",
                "right_genome",
                "coverage_similarity",
                "sampling_class",
            ]
        )
        for (left, right), (label, left_genome, right_genome, similarity, sampling_class) in sorted(
            pairs.items()
        ):
            writer.writerow(
                [
                    left,
                    right,
                    label,
                    left_genome,
                    right_genome,
                    "." if similarity is None else f"{similarity:.8f}",
                    sampling_class,
                ]
            )

    positives = sum(label == 1 for label, *_ in pairs.values())
    negatives = len(pairs) - positives
    print(
        f"bridgebin-pairs: genomes={len(by_genome)} contigs={len(records)} pairs={len(pairs)} "
        f"positive={positives} negative={negatives} coverage={'yes' if coverage else 'no'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
