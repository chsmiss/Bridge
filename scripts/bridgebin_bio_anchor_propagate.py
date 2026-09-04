#!/usr/bin/env python3
"""Propagate locally validated Biological Brain anchor groups to focused bin members.

This is the second stage of candidate-gated BridgeBin refinement. Stage one identifies
bins whose DNA anchor matrix and coverage evidence agree on a genome boundary. This stage
scores members only against those biological anchor groups and assigns by *relative*
within-sample evidence, avoiding any globally calibrated ``p_same`` threshold.

Anchor assignments are fixed by the stage-one partition. Non-anchor members use a robust
(top-k mean) score to each anchor group. Members without adequate support or with a small
between-group margin remain in the original bin unless ``--ambiguous-unbinned`` is set.
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
    p.add_argument("--assignments", type=Path, required=True)
    p.add_argument("--stage1-report", type=Path, required=True)
    p.add_argument("--focus-bins", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--top-k", type=int, default=2)
    p.add_argument("--min-group-links", type=int, default=1)
    p.add_argument("--min-margin", type=float, default=0.01)
    p.add_argument(
        "--ambiguous-unbinned",
        action="store_true",
        help="unbin ambiguous focused members instead of retaining their original bin",
    )
    return p.parse_args(argv)


def key(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def read_focus(path: Path) -> set[str]:
    result: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "bin" not in reader.fieldnames:
            raise ValueError(f"{path}: focus table requires bin column")
        for row in reader:
            value = (row.get("bin") or "").strip()
            if value:
                result.add(value)
    return result


def read_assignments(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "contig" not in reader.fieldnames or "bin" not in reader.fieldnames:
            raise ValueError(f"{path}: assignments require contig and bin columns")
        fields = list(reader.fieldnames)
        rows = [dict(row) for row in reader]
    return fields, rows


def read_scores(path: Path) -> Dict[Tuple[str, str], float]:
    result: Dict[Tuple[str, str], float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            left = (row.get("left") or row.get("source") or "").strip()
            right = (row.get("right") or row.get("target") or "").strip()
            raw = row.get("p_same") or row.get("score") or ""
            try:
                value = float(raw)
            except (TypeError, ValueError):
                continue
            if left and right and left != right and math.isfinite(value):
                result[key(left, right)] = value
    return result


def robust_group_score(values: List[float], top_k: int) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values, reverse=True)
    keep = ordered[: max(1, min(top_k, len(ordered)))]
    return sum(keep) / len(keep)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.top_k < 1 or args.min_group_links < 1:
        raise SystemExit("--top-k and --min-group-links must be >=1")

    focus = read_focus(args.focus_bins)
    stage1 = json.loads(args.stage1_report.read_text(encoding="utf-8"))
    groups: Dict[str, Tuple[List[str], List[str]]] = {}
    for entry in stage1.get("bins", []):
        name = str(entry.get("bin", ""))
        if name not in focus:
            continue
        left = [str(x) for x in entry.get("left_anchors") or []]
        right = [str(x) for x in entry.get("right_anchors") or []]
        if left and right:
            groups[name] = (left, right)

    fields, rows = read_assignments(args.assignments)
    pair_scores = read_scores(args.scores)
    members: Dict[str, List[str]] = defaultdict(list)
    row_by_contig: Dict[str, Dict[str, str]] = {}
    for row in rows:
        contig = (row.get("contig") or "").strip()
        bin_name = (row.get("bin") or "").strip()
        if contig:
            row_by_contig[contig] = row
        if contig and bin_name in groups:
            members[bin_name].append(contig)

    report_bins = []
    total_assigned = total_ambiguous = total_missing = 0
    for bin_name in sorted(groups):
        left_anchors, right_anchors = groups[bin_name]
        left_set = set(left_anchors)
        right_set = set(right_anchors)
        out_left = f"{bin_name}__bioA"
        out_right = f"{bin_name}__bioB"
        assigned_left = assigned_right = ambiguous = missing = 0
        member_rows = []

        for contig in members.get(bin_name, []):
            if contig in left_set:
                row_by_contig[contig]["bin"] = out_left
                assigned_left += 1
                member_rows.append({"contig": contig, "decision": "anchor_A", "score_a": None, "score_b": None, "margin": None})
                continue
            if contig in right_set:
                row_by_contig[contig]["bin"] = out_right
                assigned_right += 1
                member_rows.append({"contig": contig, "decision": "anchor_B", "score_a": None, "score_b": None, "margin": None})
                continue

            scores_a = [pair_scores[key(contig, anchor)] for anchor in left_anchors if key(contig, anchor) in pair_scores]
            scores_b = [pair_scores[key(contig, anchor)] for anchor in right_anchors if key(contig, anchor) in pair_scores]
            score_a = robust_group_score(scores_a, args.top_k) if len(scores_a) >= args.min_group_links else None
            score_b = robust_group_score(scores_b, args.top_k) if len(scores_b) >= args.min_group_links else None
            if score_a is None or score_b is None:
                missing += 1
                total_missing += 1
                if args.ambiguous_unbinned:
                    row_by_contig[contig]["bin"] = "unbinned"
                member_rows.append({"contig": contig, "decision": "missing", "score_a": score_a, "score_b": score_b, "margin": None})
                continue

            margin = abs(score_a - score_b)
            if margin < args.min_margin:
                ambiguous += 1
                total_ambiguous += 1
                if args.ambiguous_unbinned:
                    row_by_contig[contig]["bin"] = "unbinned"
                member_rows.append({"contig": contig, "decision": "ambiguous", "score_a": score_a, "score_b": score_b, "margin": margin})
            elif score_a > score_b:
                row_by_contig[contig]["bin"] = out_left
                assigned_left += 1
                total_assigned += 1
                member_rows.append({"contig": contig, "decision": "A", "score_a": score_a, "score_b": score_b, "margin": margin})
            else:
                row_by_contig[contig]["bin"] = out_right
                assigned_right += 1
                total_assigned += 1
                member_rows.append({"contig": contig, "decision": "B", "score_a": score_a, "score_b": score_b, "margin": margin})

        report_bins.append(
            {
                "bin": bin_name,
                "members": len(members.get(bin_name, [])),
                "anchors_a": left_anchors,
                "anchors_b": right_anchors,
                "assigned_a": assigned_left,
                "assigned_b": assigned_right,
                "ambiguous": ambiguous,
                "missing": missing,
                "member_decisions": member_rows,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    result = {
        "focus_bins": sorted(focus),
        "bins_propagated": len(report_bins),
        "top_k": args.top_k,
        "min_group_links": args.min_group_links,
        "min_margin": args.min_margin,
        "ambiguous_unbinned": args.ambiguous_unbinned,
        "assigned_nonanchors": total_assigned,
        "ambiguous_nonanchors": total_ambiguous,
        "missing_nonanchors": total_missing,
        "bins": report_bins,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"bridgebin-bio-propagate: bins={len(report_bins)} assigned={total_assigned} "
        f"ambiguous={total_ambiguous} missing={total_missing}"
    )
    for entry in report_bins:
        print(
            f"  {entry['bin']} members={entry['members']} "
            f"A={entry['assigned_a']} B={entry['assigned_b']} "
            f"ambiguous={entry['ambiguous']} missing={entry['missing']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
