#!/usr/bin/env python3
"""Detect BridgeBin bins that deserve second-pass Biological Brain expansion.

This is a production, truth-free gate between the cheap anchor probe and expensive
member-to-anchor scoring.  It consumes current assignments plus calibrated pair scores
and flags bins with enough within-bin hard-negative evidence.

The detector is deliberately permissive: false positives only cost extra model inference.
The downstream signed-cut optimizer remains responsible for deciding whether a coherent
split actually exists, so this stage never changes bin membership itself.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--pair-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="TSV of conflicted bins")
    parser.add_argument("--summary", type=Path, help="optional JSON summary")
    parser.add_argument("--split-max-same", type=float, required=True)
    parser.add_argument("--min-confidence", type=float, default=0.0)
    parser.add_argument("--min-hard-pairs", type=int, default=2)
    parser.add_argument("--min-distinct-contigs", type=int, default=2)
    return parser.parse_args(argv)


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        yield from reader


def first(row: Dict[str, str], names: Sequence[str], default: str = "") -> str:
    for name in names:
        raw = row.get(name)
        if raw is not None and raw.strip() and raw.strip() not in {".", "NA", "na"}:
            return raw.strip()
    return default


def read_assignments(path: Path) -> Dict[str, Optional[str]]:
    result: Dict[str, Optional[str]] = {}
    for row in rows(path):
        contig = first(row, ("contig", "contig_id", "sequence"))
        raw_bin = first(row, ("bin", "bin_id", "cluster"))
        if not contig:
            continue
        result[contig] = None if raw_bin in {"", ".", "NA", "unbinned"} else raw_bin
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.split_max_same <= 1.0:
        raise ValueError("--split-max-same must be in [0,1]")
    if not 0.0 <= args.min_confidence <= 1.0:
        raise ValueError("--min-confidence must be in [0,1]")
    if args.min_hard_pairs < 1 or args.min_distinct_contigs < 2:
        raise ValueError("support limits must be positive and need at least two contigs")

    assignments = read_assignments(args.assignments)
    evidence = defaultdict(lambda: {"pairs": 0, "contigs": set(), "scores": []})
    scored_pairs = 0
    same_bin_pairs = 0

    for row in rows(args.pair_scores):
        left = first(row, ("left", "source", "contig_a", "contig1"))
        right = first(row, ("right", "target", "contig_b", "contig2"))
        if not left or not right or left == right:
            continue
        raw_same = first(row, ("p_same", "same_probability", "same_genome", "probability", "score"))
        if not raw_same:
            continue
        p_same = float(raw_same)
        if not math.isfinite(p_same) or not 0.0 <= p_same <= 1.0:
            raise ValueError(f"invalid p_same for {left!r},{right!r}: {raw_same!r}")
        raw_conf = first(row, ("confidence", "model_confidence", "pair_confidence"))
        confidence = 1.0 if not raw_conf else float(raw_conf)
        if not math.isfinite(confidence) or not 0.0 <= confidence <= 1.0:
            raise ValueError(f"invalid confidence for {left!r},{right!r}: {raw_conf!r}")
        scored_pairs += 1
        left_bin = assignments.get(left)
        right_bin = assignments.get(right)
        if left_bin is None or left_bin != right_bin:
            continue
        same_bin_pairs += 1
        if confidence < args.min_confidence or p_same > args.split_max_same:
            continue
        record = evidence[left_bin]
        record["pairs"] += 1
        record["contigs"].update((left, right))
        record["scores"].append(p_same)

    selected = []
    for bin_name, record in sorted(evidence.items()):
        pair_count = int(record["pairs"])
        contigs = record["contigs"]
        scores = record["scores"]
        if pair_count < args.min_hard_pairs or len(contigs) < args.min_distinct_contigs:
            continue
        selected.append(
            {
                "bin": bin_name,
                "hard_pairs": pair_count,
                "distinct_contigs": len(contigs),
                "min_p_same": min(scores),
                "mean_p_same": sum(scores) / len(scores),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("bin", "hard_pairs", "distinct_contigs", "min_p_same", "mean_p_same"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in selected:
            writer.writerow(
                {
                    **row,
                    "min_p_same": f"{row['min_p_same']:.8f}",
                    "mean_p_same": f"{row['mean_p_same']:.8f}",
                }
            )

    summary = {
        "scored_pairs": scored_pairs,
        "same_bin_pairs": same_bin_pairs,
        "bins_with_any_hard_pair": len(evidence),
        "selected_conflicted_bins": len(selected),
        "selected_bins": [row["bin"] for row in selected],
        "split_max_same": args.split_max_same,
        "min_confidence": args.min_confidence,
        "min_hard_pairs": args.min_hard_pairs,
        "min_distinct_contigs": args.min_distinct_contigs,
    }
    if args.summary is not None:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
