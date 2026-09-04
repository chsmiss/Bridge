#!/usr/bin/env python3
"""Select bins for expensive Biological Brain expansion using aggregate conflict consensus.

A calibrated pair head can transfer useful ranking while its absolute probabilities shift
between datasets.  Treating every low pair as a cannot-link therefore over-splits pure
bins.  This gate promotes a bin only when low same-genome probabilities are a *bin-level*
phenomenon across the probe anchor graph.

The procedure is truth-free: it uses current assignments, candidate classes, pair scores,
and the model's reference split threshold.  A single low pair can never select a bin.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assignments", type=Path, required=True)
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--split-threshold", type=float, required=True)
    p.add_argument("--min-anchor-pairs", type=int, default=15)
    p.add_argument("--min-conflict-fraction", type=float, default=0.50)
    p.add_argument("--min-bin-bp", type=int, default=250000)
    p.add_argument("--min-score-confidence", type=float, default=0.0)
    return p.parse_args(argv)


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        yield from reader


def key(left: str, right: str) -> Tuple[str, str]:
    return (left, right) if left <= right else (right, left)


def read_assignments(path: Path) -> Tuple[Dict[str, Optional[str]], Dict[str, int], Dict[str, int]]:
    assignment: Dict[str, Optional[str]] = {}
    length: Dict[str, int] = {}
    bp: Dict[str, int] = defaultdict(int)
    for row in rows(path):
        contig = (row.get("contig") or "").strip()
        raw_bin = (row.get("bin") or "").strip()
        if not contig:
            continue
        bin_id = None if raw_bin in {"", ".", "NA", "na", "unbinned"} else raw_bin
        try:
            contig_length = int(float((row.get("length") or "0").strip() or "0"))
        except ValueError as error:
            raise ValueError(f"{path}: invalid length for {contig!r}") from error
        assignment[contig] = bin_id
        length[contig] = max(0, contig_length)
        if bin_id is not None:
            bp[bin_id] += max(0, contig_length)
    return assignment, length, dict(bp)


def read_scores(path: Path) -> Dict[Tuple[str, str], Tuple[float, float]]:
    out: Dict[Tuple[str, str], Tuple[float, float]] = {}
    for row in rows(path):
        left = (row.get("left") or row.get("source") or "").strip()
        right = (row.get("right") or row.get("target") or "").strip()
        if not left or not right or left == right:
            continue
        p = float(row.get("p_same") or row.get("score") or "nan")
        confidence = float(row.get("confidence") or "1")
        if math.isfinite(p) and 0.0 <= p <= 1.0 and math.isfinite(confidence):
            out[key(left, right)] = (p, max(0.0, min(1.0, confidence)))
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.split_threshold <= 1.0:
        raise SystemExit("--split-threshold must be in [0,1]")
    if args.min_anchor_pairs < 1 or args.min_bin_bp < 0:
        raise SystemExit("pair/bp minimums must be non-negative")
    if not 0.0 < args.min_conflict_fraction <= 1.0:
        raise SystemExit("--min-conflict-fraction must be in (0,1]")
    if not 0.0 <= args.min_score_confidence <= 1.0:
        raise SystemExit("--min-score-confidence must be in [0,1]")

    assignment, _length, bin_bp = read_assignments(args.assignments)
    score = read_scores(args.scores)
    by_bin: Dict[str, List[float]] = defaultdict(list)

    for row in rows(args.candidates):
        classes = {
            piece.strip()
            for piece in (row.get("candidate_class") or "").split(";")
            if piece.strip()
        }
        if "within_bin_anchor" not in classes:
            continue
        left = (row.get("left") or "").strip()
        right = (row.get("right") or "").strip()
        left_bin = assignment.get(left)
        right_bin = assignment.get(right)
        if left_bin is None or left_bin != right_bin:
            continue
        value = score.get(key(left, right))
        if value is None or value[1] < args.min_score_confidence:
            continue
        by_bin[left_bin].append(value[0])

    selected = []
    diagnostics = []
    for bin_id in sorted(by_bin):
        values = by_bin[bin_id]
        hard = sum(value <= args.split_threshold for value in values)
        fraction = hard / max(1, len(values))
        median = statistics.median(values)
        bp = bin_bp.get(bin_id, 0)
        focus = (
            len(values) >= args.min_anchor_pairs
            and bp >= args.min_bin_bp
            and fraction >= args.min_conflict_fraction
        )
        diagnostics.append((bin_id, bp, len(values), hard, fraction, median, focus))
        if focus:
            selected.append(bin_id)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "bin",
                "bp",
                "anchor_pairs",
                "conflict_pairs",
                "conflict_fraction",
                "median_p_same",
                "focus",
            ]
        )
        for bin_id, bp, pairs, hard, fraction, median, focus in diagnostics:
            if focus:
                writer.writerow(
                    [bin_id, bp, pairs, hard, f"{fraction:.8f}", f"{median:.8f}", 1]
                )

    print(
        f"bridgebin-focus-conflicts: probed_bins={len(diagnostics)} focus_bins={len(selected)} "
        f"split_threshold={args.split_threshold:.6g} min_conflict_fraction={args.min_conflict_fraction:.3f} "
        f"min_bin_bp={args.min_bin_bp}"
    )
    if selected:
        print("bridgebin-focus-conflicts: selected=" + ",".join(selected))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
