#!/usr/bin/env python3
"""Truth-free local calibration of Biological Brain anchor scores.

The absolute probability emitted by a reference-trained same-genome head may shift across
communities. This probe therefore ignores global probability thresholds. For each current
bin it uses only the complete-ish ``within_bin_anchor`` score graph and asks whether the
anchor matrix contains a statistically coherent two-block structure: within-block scores
must exceed cross-block scores by more than expected after edge-score permutation.

A DNA-only block can reflect chromosome-region or strain-like representation structure
inside one genome, so a split is accepted only when sample-specific coverage supports the
same anchor partition. Biological evidence proposes the identity boundary; independent
coverage evidence validates that it behaves like a genome boundary in this sample. The
remaining contigs are propagated to the accepted anchor groups with the cheap
coverage/composition/GC similarity already used by BridgeBin candidate mining.

This script is intended as an experimental pre-refiner. It never uses benchmark truth.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import bridgebin_candidate_pairs as cheap


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs", type=Path, required=True)
    p.add_argument("--coverage", type=Path)
    p.add_argument("--assignments", type=Path, required=True)
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--min-length", type=int, default=1500)
    p.add_argument("--min-anchors", type=int, default=6)
    p.add_argument("--min-side-anchors", type=int, default=2)
    p.add_argument("--min-edge-density", type=float, default=0.80)
    p.add_argument("--min-gap", type=float, default=0.035)
    p.add_argument("--max-permutation-p", type=float, default=0.05)
    p.add_argument(
        "--min-coverage-gap",
        type=float,
        default=0.0,
        help=(
            "require mean within-group coverage similarity minus cross-group coverage "
            "similarity to be at least this value; 0 disables the extra consensus gate"
        ),
    )
    p.add_argument("--permutations", type=int, default=96)
    p.add_argument("--member-margin", type=float, default=0.0)
    p.add_argument("--seed", type=int, default=43)
    return p.parse_args(argv)


def pair_key(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def read_assignment_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]], Dict[str, Optional[str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "contig" not in reader.fieldnames or "bin" not in reader.fieldnames:
            raise ValueError(f"{path}: assignments require contig and bin columns")
        rows = [dict(row) for row in reader]
        fields = list(reader.fieldnames)
    mapping: Dict[str, Optional[str]] = {}
    for row in rows:
        contig = (row.get("contig") or "").strip()
        raw = (row.get("bin") or "").strip()
        if contig:
            mapping[contig] = None if raw in {"", ".", "NA", "unbinned"} else raw
    return fields, rows, mapping


def read_anchor_pairs(path: Path) -> Dict[Tuple[str, str], str]:
    result: Dict[Tuple[str, str], str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            classes = (row.get("candidate_class") or "").split(";")
            if "within_bin_anchor" not in classes:
                continue
            left = (row.get("left") or "").strip()
            right = (row.get("right") or "").strip()
            if left and right and left != right:
                result[pair_key(left, right)] = "within_bin_anchor"
    return result


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
                result[pair_key(left, right)] = value
    return result


def mean(values: Iterable[float]) -> float:
    data = list(values)
    return sum(data) / len(data) if data else float("nan")


def evaluate_partition(
    anchors: List[str],
    edges: Dict[Tuple[str, str], float],
    mask: int,
) -> Optional[Tuple[float, float, float, List[str], List[str]]]:
    left = [anchors[0]]
    right: List[str] = []
    for idx in range(1, len(anchors)):
        (right if (mask >> (idx - 1)) & 1 else left).append(anchors[idx])
    if len(left) < 2 or len(right) < 2:
        return None
    within: List[float] = []
    cross: List[float] = []
    side = {name: 0 for name in left}
    side.update({name: 1 for name in right})
    for (a, b), value in edges.items():
        if a not in side or b not in side:
            continue
        (within if side[a] == side[b] else cross).append(value)
    if not within or not cross:
        return None
    within_mean = mean(within)
    cross_mean = mean(cross)
    return within_mean - cross_mean, within_mean, cross_mean, left, right


def best_partition(
    anchors: List[str],
    edges: Dict[Tuple[str, str], float],
    min_side: int,
) -> Optional[Tuple[float, float, float, List[str], List[str]]]:
    best = None
    for mask in range(1 << (len(anchors) - 1)):
        candidate = evaluate_partition(anchors, edges, mask)
        if candidate is None or min(len(candidate[3]), len(candidate[4])) < min_side:
            continue
        if best is None or candidate[0] > best[0] + 1e-12:
            best = candidate
    return best


def permuted_p_value(
    anchors: List[str],
    edges: Dict[Tuple[str, str], float],
    observed_gap: float,
    min_side: int,
    permutations: int,
    rng: random.Random,
) -> float:
    keys = list(edges)
    values = [edges[key] for key in keys]
    exceed = 0
    for _ in range(max(0, permutations)):
        shuffled = list(values)
        rng.shuffle(shuffled)
        null_edges = dict(zip(keys, shuffled))
        null_best = best_partition(anchors, null_edges, min_side)
        if null_best is not None and null_best[0] >= observed_gap - 1e-12:
            exceed += 1
    return (exceed + 1.0) / (max(0, permutations) + 1.0)


def coverage_partition_support(
    left: List[cheap.Feature], right: List[cheap.Feature]
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    within: List[float] = []
    cross: List[float] = []
    for group in (left, right):
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                coverage, _composition, _gc = cheap.component_sim(group[i], group[j])
                if coverage is not None:
                    within.append(coverage)
    for a in left:
        for b in right:
            coverage, _composition, _gc = cheap.component_sim(a, b)
            if coverage is not None:
                cross.append(coverage)
    if not within or not cross:
        return None, None, None
    within_mean = mean(within)
    cross_mean = mean(cross)
    return within_mean - cross_mean, within_mean, cross_mean


def group_similarity(member: cheap.Feature, anchors: List[cheap.Feature]) -> float:
    scores = sorted((cheap.cheap_similarity(member, anchor) for anchor in anchors), reverse=True)
    if not scores:
        return 0.0
    keep = scores[: min(3, len(scores))]
    return sum(keep) / len(keep)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    fields, rows, assignments = read_assignment_rows(args.assignments)
    anchor_pairs = read_anchor_pairs(args.candidates)
    raw_scores = read_scores(args.scores)
    coverage = cheap.read_coverage(args.coverage)
    features = cheap.build_features(args.contigs, coverage, args.min_length)

    bins: Dict[str, List[str]] = defaultdict(list)
    for contig, bin_name in assignments.items():
        if bin_name is not None and contig in features:
            bins[bin_name].append(contig)

    bin_edges: Dict[str, Dict[Tuple[str, str], float]] = defaultdict(dict)
    bin_anchors: Dict[str, set[str]] = defaultdict(set)
    for key in anchor_pairs:
        left, right = key
        left_bin = assignments.get(left)
        right_bin = assignments.get(right)
        if left_bin is None or left_bin != right_bin or key not in raw_scores:
            continue
        bin_edges[left_bin][key] = raw_scores[key]
        bin_anchors[left_bin].update(key)

    rng = random.Random(args.seed)
    report_bins = []
    rewrites: Dict[str, Optional[str]] = {}
    split_count = 0
    ambiguous = 0

    for bin_name in sorted(bins):
        anchors = sorted(bin_anchors.get(bin_name, set()))
        edges = bin_edges.get(bin_name, {})
        possible = len(anchors) * (len(anchors) - 1) // 2
        density = len(edges) / possible if possible else 0.0
        entry = {
            "bin": bin_name,
            "anchors": len(anchors),
            "edges": len(edges),
            "edge_density": density,
            "split": False,
        }
        if len(anchors) < args.min_anchors or density < args.min_edge_density:
            report_bins.append(entry)
            continue
        best = best_partition(anchors, edges, args.min_side_anchors)
        if best is None:
            report_bins.append(entry)
            continue
        gap, within_mean, cross_mean, left_anchor_ids, right_anchor_ids = best
        p_value = permuted_p_value(
            anchors, edges, gap, args.min_side_anchors, args.permutations, rng
        )
        entry.update(
            {
                "gap": gap,
                "within_mean": within_mean,
                "cross_mean": cross_mean,
                "permutation_p": p_value,
                "left_anchors": left_anchor_ids,
                "right_anchors": right_anchor_ids,
            }
        )
        if gap < args.min_gap or p_value > args.max_permutation_p:
            report_bins.append(entry)
            continue

        left_features = [features[name] for name in left_anchor_ids if name in features]
        right_features = [features[name] for name in right_anchor_ids if name in features]
        if len(left_features) < args.min_side_anchors or len(right_features) < args.min_side_anchors:
            report_bins.append(entry)
            continue

        coverage_gap, coverage_within, coverage_cross = coverage_partition_support(
            left_features, right_features
        )
        entry.update(
            {
                "coverage_gap": coverage_gap,
                "coverage_within_mean": coverage_within,
                "coverage_cross_mean": coverage_cross,
            }
        )
        if args.min_coverage_gap > 0.0 and (
            coverage_gap is None or coverage_gap < args.min_coverage_gap
        ):
            report_bins.append(entry)
            continue

        left_bin = f"{bin_name}__bioA"
        right_bin = f"{bin_name}__bioB"
        assigned_left = assigned_right = unassigned = 0
        for contig in bins[bin_name]:
            feature = features[contig]
            left_score = group_similarity(feature, left_features)
            right_score = group_similarity(feature, right_features)
            margin = abs(left_score - right_score)
            if margin < args.member_margin:
                rewrites[contig] = None
                unassigned += 1
                ambiguous += 1
            elif left_score >= right_score:
                rewrites[contig] = left_bin
                assigned_left += 1
            else:
                rewrites[contig] = right_bin
                assigned_right += 1
        entry.update(
            {
                "split": True,
                "left_bin": left_bin,
                "right_bin": right_bin,
                "assigned_left": assigned_left,
                "assigned_right": assigned_right,
                "unassigned_margin": unassigned,
            }
        )
        split_count += 1
        report_bins.append(entry)

    for row in rows:
        contig = (row.get("contig") or "").strip()
        if contig not in rewrites:
            continue
        new_bin = rewrites[contig]
        row["bin"] = "unbinned" if new_bin is None else new_bin

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "bins_considered": len(bins),
        "bins_split": split_count,
        "ambiguous_unbinned": ambiguous,
        "min_gap": args.min_gap,
        "max_permutation_p": args.max_permutation_p,
        "min_coverage_gap": args.min_coverage_gap,
        "permutations": args.permutations,
        "member_margin": args.member_margin,
        "bins": report_bins,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"bridgebin-local-anchor: bins={len(bins)} split={split_count} "
        f"ambiguous_unbinned={ambiguous}"
    )
    for entry in report_bins:
        if entry.get("split"):
            coverage_text = (
                "." if entry.get("coverage_gap") is None else f"{entry['coverage_gap']:.6f}"
            )
            print(
                "  split "
                f"{entry['bin']} anchors={entry['anchors']} gap={entry['gap']:.6f} "
                f"p={entry['permutation_p']:.6f} coverage_gap={coverage_text} "
                f"sizes={entry['assigned_left']}/{entry['assigned_right']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
