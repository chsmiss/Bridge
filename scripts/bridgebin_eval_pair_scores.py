#!/usr/bin/env python3
"""Evaluate same-genome pair probabilities against a labeled pair table.

This is benchmark/training diagnostics only.  It reports threshold-free AUC, confusion
statistics at supplied join/split thresholds, and optional focus-pair distributions such
as Escherichia-vs-Salmonella.  Production scoring remains truth-free.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
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
        lo = int(math.floor(pos)); hi = int(math.ceil(pos))
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


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.split_threshold < args.join_threshold <= 1.0:
        raise SystemExit("require 0 <= split-threshold < join-threshold <= 1")

    truth: Dict[Tuple[str, str], Tuple[int, str, str]] = {}
    with args.pairs.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            raw = (row.get("label") or row.get("same_genome") or "").strip()
            if not left or not right or raw not in {"0", "1", "0.0", "1.0"}:
                continue
            truth[key(left, right)] = (
                int(float(raw)),
                (row.get("left_genome") or row.get("left_group") or "").strip(),
                (row.get("right_genome") or row.get("right_group") or "").strip(),
            )

    scored: List[Tuple[float, int, str, str]] = []
    with args.scores.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            if key(left, right) not in truth:
                continue
            p = float(row.get("p_same") or row.get("score") or "nan")
            if not math.isfinite(p):
                continue
            label, lg, rg = truth[key(left, right)]
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
                    p for p, _y, lg, rg in scored
                    if {lg, rg} == {left, right}
                ]
                focus[f"{left} vs {right}"] = quantiles(values)

    result = {
        "scored_pairs": len(scored),
        "positive_pairs": sum(y for _, y, _, _ in scored),
        "negative_pairs": sum(1 - y for _, y, _, _ in scored),
        "auc": auc([(p, y) for p, y, _, _ in scored]),
        "same_genome_scores": quantiles(positives),
        "different_genome_scores": quantiles(negatives),
        "join_threshold": args.join_threshold,
        "join_pairs": len(join),
        "join_precision": None if not join else sum(y for _, y, _, _ in join) / len(join),
        "join_recall": None if not positives else sum(y for _, y, _, _ in join) / len(positives),
        "split_threshold": args.split_threshold,
        "split_pairs": len(split),
        "split_precision": None if not split else sum(1 - y for _, y, _, _ in split) / len(split),
        "split_recall": None if not negatives else sum(1 - y for _, y, _, _ in split) / len(negatives),
        "ambiguous_pairs": len(ambiguous),
        "focus": focus,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
