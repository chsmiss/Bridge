#!/usr/bin/env python3
"""Benchmark-only evaluation of BridgeBin Biological Brain candidate mining.

Production candidate mining is truth-free. This script joins its emitted pairs to known
benchmark genomes and asks whether the expensive pair model is being shown the decisions
that matter:

- within-bin hard-negative exposure: mixed current bins need cross-genome candidate pairs;
- cross-bin recovery exposure: fragmented genomes need same-genome candidate pairs across bins;
- residual rescue exposure: an unbinned contig needs at least one same-genome candidate;
- pair budget / class composition.

Truth is used only here, never by ``bridgebin_candidate_pairs.py``.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def read_truth(path: Path) -> Dict[str, Tuple[str, int]]:
    result: Dict[str, Tuple[str, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            genome = (row.get("genome") or row.get("species") or "").strip()
            if not contig or not genome:
                continue
            eligible = (row.get("eligible") or "1").strip().lower() not in {
                "0",
                "false",
                "no",
            }
            if not eligible:
                continue
            length = int(float(row.get("length") or 0))
            result[contig] = (genome, length)
    return result


def read_assignments(path: Path) -> Dict[str, Optional[str]]:
    result: Dict[str, Optional[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "contig" not in reader.fieldnames or "bin" not in reader.fieldnames:
            raise ValueError(f"{path}: assignments need contig and bin columns")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            raw_bin = (row.get("bin") or "").strip()
            if not contig:
                continue
            result[contig] = None if raw_bin in {"", ".", "NA", "unbinned"} else raw_bin
    return result


def read_candidates(path: Path) -> List[Tuple[str, str, Set[str]]]:
    result = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            if not left or not right or left == right:
                continue
            classes = {
                value.strip()
                for value in (row.get("candidate_class") or "").split(";")
                if value.strip()
            }
            result.append((left, right, classes))
    return result


def safe_fraction(numerator: int, denominator: int) -> Optional[float]:
    return None if denominator == 0 else numerator / denominator


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    truth = read_truth(args.truth)
    assignments = read_assignments(args.assignments)
    candidates = read_candidates(args.candidates)

    by_bin: Dict[str, List[str]] = defaultdict(list)
    residuals = []
    genome_bins: Dict[str, Set[str]] = defaultdict(set)
    for contig, (genome, _length) in truth.items():
        bin_name = assignments.get(contig)
        if bin_name is None:
            residuals.append(contig)
        else:
            by_bin[bin_name].append(contig)
            genome_bins[genome].add(bin_name)

    mixed_bins = {
        bin_name
        for bin_name, members in by_bin.items()
        if len({truth[contig][0] for contig in members}) > 1
    }
    fragmented_genomes = {
        genome for genome, bins in genome_bins.items() if len(bins) > 1
    }

    mixed_bins_exposed: Set[str] = set()
    fragmented_genomes_exposed: Set[str] = set()
    residuals_exposed: Set[str] = set()
    same_pairs = 0
    different_pairs = 0
    class_counts: Dict[str, int] = defaultdict(int)
    same_class_counts: Dict[str, int] = defaultdict(int)
    different_class_counts: Dict[str, int] = defaultdict(int)

    for left, right, classes in candidates:
        if left not in truth or right not in truth:
            continue
        left_genome = truth[left][0]
        right_genome = truth[right][0]
        same = left_genome == right_genome
        left_bin = assignments.get(left)
        right_bin = assignments.get(right)
        if same:
            same_pairs += 1
        else:
            different_pairs += 1
        for candidate_class in classes:
            class_counts[candidate_class] += 1
            (same_class_counts if same else different_class_counts)[candidate_class] += 1

        if (
            not same
            and left_bin is not None
            and left_bin == right_bin
            and left_bin in mixed_bins
        ):
            mixed_bins_exposed.add(left_bin)
        if (
            same
            and left_bin is not None
            and right_bin is not None
            and left_bin != right_bin
            and left_genome in fragmented_genomes
        ):
            fragmented_genomes_exposed.add(left_genome)
        if same:
            if left_bin is None and right_bin is not None:
                residuals_exposed.add(left)
            if right_bin is None and left_bin is not None:
                residuals_exposed.add(right)

    metrics = {
        "candidate_pairs": len(candidates),
        "truth_scored_pairs": same_pairs + different_pairs,
        "same_genome_pairs": same_pairs,
        "different_genome_pairs": different_pairs,
        "class_counts": dict(sorted(class_counts.items())),
        "same_genome_class_counts": dict(sorted(same_class_counts.items())),
        "different_genome_class_counts": dict(sorted(different_class_counts.items())),
        "mixed_bins": len(mixed_bins),
        "mixed_bins_with_cross_genome_candidate": len(mixed_bins_exposed),
        "mixed_bin_exposure_recall": safe_fraction(len(mixed_bins_exposed), len(mixed_bins)),
        "fragmented_genomes": len(fragmented_genomes),
        "fragmented_genomes_with_cross_bin_same_candidate": len(fragmented_genomes_exposed),
        "fragmented_genome_exposure_recall": safe_fraction(
            len(fragmented_genomes_exposed), len(fragmented_genomes)
        ),
        "eligible_residual_contigs": len(residuals),
        "residuals_with_same_genome_candidate": len(residuals_exposed),
        "residual_exposure_recall": safe_fraction(len(residuals_exposed), len(residuals)),
        "mixed_bins_exposed": sorted(mixed_bins_exposed),
        "fragmented_genomes_exposed": sorted(fragmented_genomes_exposed),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(metrics, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
