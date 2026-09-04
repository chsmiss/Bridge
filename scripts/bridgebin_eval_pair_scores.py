#!/usr/bin/env python3
"""Evaluate same-genome pair probabilities against benchmark truth.

Two diagnostic modes are supported:

1. a pair table that already contains ``label``/``same_genome`` plus optional genome
   columns (training/validation diagnostics);
2. an unlabeled production-style candidate table plus ``--truth`` mapping contigs to
   benchmark genomes (end-to-end benchmark diagnostics).

This script is never used for production decisions. It reports threshold-free AUC,
confusion statistics at supplied join/split thresholds, and optional focus-pair
score distributions such as Escherichia-vs-Salmonella.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--truth", type=Path, help="optional contig->genome benchmark truth table")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--join-threshold", type=float, default=0.88)
    p.add_argument("--split-threshold", type=float, default=0.12)
    p.add_argument("--focus", default="")
    return p.parse_args(argv)


def key(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a < b else (b, a)


def quantiles(values: List[float]) -> Dict[str, Optional[float]]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "q05": None, "q95": None}
    x = sorted(values)

    def q(f: float) -> float:
        pos = f * (len(x) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return x[lo]
        w = pos - lo
        return x[lo] * (1.0 - w) + x[hi] * w

    return {
        "n": len(x),
        "mean": sum(x) / len(x),
        "median": q(0.5),
        "q05": q(0.05),
        "q95": q(0.95),
    }


def auc(rows: List[Tuple[float, int]]) -> Optional[float]:
    pos = [p for p, y in rows if y == 1]
    neg = [p for p, y in rows if y == 0]
    if not pos or not neg:
        return None
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else 0.5 if p == n else 0.0
    return wins / (len(pos) * len(neg))


def read_contig_truth(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    result: Dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            genome = (row.get("genome") or row.get("species") or row.get("group") or "").strip()
            eligible = (row.get("eligible") or "1").strip().lower() not in {
                "0",
                "false",
                "no",
            }
            if contig and genome and eligible:
                result[contig] = genome
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.split_threshold < args.join_threshold <= 1.0:
        raise SystemExit("require 0 <= split-threshold < join-threshold <= 1")

    contig_truth = read_contig_truth(args.truth)
    truth: Dict[Tuple[str, str], Tuple[int, str, str]] = {}
    with args.pairs.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            if not left or not right or left == right:
                continue
            raw = (row.get("label") or row.get("same_genome") or "").strip()
            if raw in {"0", "1", "0.0", "1.0"}:
                truth[key(left, right)] = (
                    int(float(raw)),
                    (row.get("left_genome") or row.get("left_group") or "").strip(),
                    (row.get("right_genome") or row.get("right_group") or "").strip(),
                )
                continue
            if left in contig_truth and right in contig_truth:
                lg = contig_truth[left]
                rg = contig_truth[right]
                truth[key(left, right)] = (int(lg == rg), lg, rg)

    scored: List[Tuple[float, int, str, str]] = []
    with args.scores.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            pair_key = key(left, right)
            if pair_key not in truth:
                continue
            p = float(row.get("p_same") or row.get("score") or "nan")
            if not math.isfinite(p):
                continue
            label, lg, rg = truth[pair_key]
            scored.append((p, label, lg, rg))
    if not scored:
        raise SystemExit("no scored labeled pairs")

    join = [row for row in scored if row[0] >= args.join_threshold]
    split = [row for row in scored if row[0] <= args.split_threshold]
    ambiguous = [row for row in scored if args.split_threshold < row[0] < args.join_threshold]
    positives = [p for p, y, _, _ in scored if y == 1]
    negatives = [p for p, y, _, _ in scored if y == 0]

    focus_groups = [x.strip() for x in args.focus.split(",") if x.strip()]
    focus: Dict[str, Dict[str, Optional[float]]] = {}
    if len(focus_groups) >= 2:
        for i, left in enumerate(focus_groups):
            for right in focus_groups[i + 1 :]:
                values = [
                    p
                    for p, _y, lg, rg in scored
                    if {lg, rg} == {left, right}
                ]
                focus[f"{left} vs {right}"] = quantiles(values)

    result = {
        "scored_pairs": len(scored),
        "positive_pairs": len(positives),
        "negative_pairs": len(negatives),
        "auc": auc([(p, y) for p, y, _, _ in scored]),
        "join_threshold": args.join_threshold,
        "split_threshold": args.split_threshold,
        "join_pairs": len(join),
        "split_pairs": len(split),
        "ambiguous_pairs": len(ambiguous),
        "join_precision": None
        if not join
        else sum(y for _, y, _, _ in join) / len(join),
        "join_recall": None
        if not positives
        else sum(y for _, y, _, _ in join) / len(positives),
        "split_precision": None
        if not split
        else sum(1 - y for _, y, _, _ in split) / len(split),
        "split_recall": None
        if not negatives
        else sum(1 - y for _, y, _, _ in split) / len(negatives),
        "same_genome_scores": quantiles(positives),
        "different_genome_scores": quantiles(negatives),
        "focus": focus,
        "truth_mode": "benchmark_truth" if args.truth is not None else "labeled_pairs",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
