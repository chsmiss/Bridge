#!/usr/bin/env python3
"""Gate Biological Brain pair evidence before the Rust v2.1 consumer.

Within-bin low probabilities are allowed to drive splitting only for bins selected by the
bin-level conflict-consensus probe.  Other within-bin scores are omitted, preventing an
OOD absolute threshold from declaring every pure bin conflicted.  Cross-bin scores are
kept only when they pass the conservative learned join threshold, so safe merge evidence
can still propagate.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Set, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assignments", type=Path, required=True)
    p.add_argument("--focus-bins", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--join-threshold", type=float, required=True)
    return p.parse_args(argv)


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        yield from reader


def read_assignments(path: Path) -> Dict[str, Optional[str]]:
    out: Dict[str, Optional[str]] = {}
    for row in rows(path):
        contig = (row.get("contig") or "").strip()
        raw_bin = (row.get("bin") or "").strip()
        if not contig:
            continue
        out[contig] = None if raw_bin in {"", ".", "NA", "na", "unbinned"} else raw_bin
    return out


def read_focus(path: Path) -> Set[str]:
    focus: Set[str] = set()
    for row in rows(path):
        value = (row.get("bin") or row.get("bin_id") or "").strip()
        if value:
            focus.add(value)
    return focus


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.join_threshold <= 1.0:
        raise SystemExit("--join-threshold must be in [0,1]")
    assignment = read_assignments(args.assignments)
    focus = read_focus(args.focus_bins)

    kept_within = kept_cross = dropped = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.scores.open("r", encoding="utf-8", newline="") as source, args.output.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{args.scores}: missing header")
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            if not left or not right:
                dropped += 1
                continue
            left_bin = assignment.get(left)
            right_bin = assignment.get(right)
            raw = row.get("p_same") or row.get("score") or "nan"
            try:
                probability = float(raw)
            except ValueError:
                dropped += 1
                continue
            if not math.isfinite(probability):
                dropped += 1
                continue

            if left_bin is not None and left_bin == right_bin:
                if left_bin in focus:
                    writer.writerow(row)
                    kept_within += 1
                else:
                    dropped += 1
                continue

            if left_bin is not None and right_bin is not None and probability >= args.join_threshold:
                writer.writerow(row)
                kept_cross += 1
            else:
                dropped += 1

    print(
        f"bridgebin-filter-pair-scores: focus_bins={len(focus)} kept_within={kept_within} "
        f"kept_cross={kept_cross} dropped={dropped}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
