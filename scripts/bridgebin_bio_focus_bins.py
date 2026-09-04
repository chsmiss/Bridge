#!/usr/bin/env python3
"""Select truth-free BridgeBin bins for second-stage Biological Brain scoring.

Consumes the local anchor partition report. A bin is selected only when the DNA anchor
matrix has a significant two-block structure and independent sample-specific coverage
supports the same partition. The output is a one-column TSV accepted by
``bridgebin_candidate_pairs.py --focus-bins``.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Optional, Sequence


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--min-dna-gap", type=float, default=0.035)
    p.add_argument("--max-permutation-p", type=float, default=0.05)
    p.add_argument("--min-coverage-gap", type=float, default=0.03)
    p.add_argument(
        "--singleton-min-coverage-gap",
        type=float,
        default=0.05,
        help="stricter coverage support required when one biological anchor group has one anchor",
    )
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    selected = []
    for entry in report.get("bins", []):
        gap = entry.get("gap")
        p_value = entry.get("permutation_p")
        coverage_gap = entry.get("coverage_gap")
        left = entry.get("left_anchors") or []
        right = entry.get("right_anchors") or []
        if gap is None or p_value is None or coverage_gap is None:
            continue
        if gap < args.min_dna_gap or p_value > args.max_permutation_p:
            continue
        required_coverage = args.min_coverage_gap
        if min(len(left), len(right)) <= 1:
            required_coverage = max(required_coverage, args.singleton_min_coverage_gap)
        if coverage_gap < required_coverage:
            continue
        selected.append(
            (
                str(entry["bin"]),
                float(gap),
                float(p_value),
                float(coverage_gap),
                len(left),
                len(right),
            )
        )

    selected.sort(key=lambda row: (-row[3], -row[1], row[0]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["bin", "dna_gap", "permutation_p", "coverage_gap", "left_anchors", "right_anchors"])
        writer.writerows(selected)

    summary = ",".join(row[0] for row in selected) if selected else "none"
    print(f"bridgebin-bio-focus: selected={len(selected)} bins={summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
