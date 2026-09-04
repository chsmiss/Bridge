#!/usr/bin/env python3
"""Target-specific few-shot metric adaptation for suspicious BridgeBin bins.

Stage-1 Biological Brain already supplies a *truth-free* two-group partition of diverse
anchors, and ``bridgebin_bio_focus_bins.py`` keeps only partitions independently supported
by sample-specific evidence.  This script treats those anchor groups as local pseudo-labels
and trains a tiny residual adapter per focused bin on top of frozen DNABERT-S contig
embeddings.  It then assigns all members by relative similarity to adapted group
prototypes.  No cross-dataset absolute probability threshold is used.

The adapter starts as the identity and is regularized toward the pretrained geometry, so
bad or weak pseudo-labels cannot arbitrarily rewrite the whole representation space.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import torch
import torch.nn.functional as F

from bridgebin_metric_transfer import ResidualMetricAdapter


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", type=Path, required=True)
    p.add_argument("--assignments", type=Path, required=True)
    p.add_argument("--stage1-report", type=Path, required=True)
    p.add_argument("--focus-bins", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--rank", type=int, default=32)
    p.add_argument("--epochs", type=int, default=250)
    p.add_argument("--learning-rate", type=float, default=5e-3)
    p.add_argument("--negative-margin", type=float, default=0.20)
    p.add_argument("--negative-weight", type=float, default=5.0)
    p.add_argument("--drift-weight", type=float, default=0.15)
    p.add_argument("--residual-scale", type=float, default=0.30)
    p.add_argument("--member-margin", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=43)
    return p.parse_args(argv)


def parse_vector(raw: str) -> List[float]:
    values = [float(piece) for piece in raw.split(",") if piece.strip()]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("invalid DNA embedding")
    return values


def read_features(path: Path) -> Dict[str, List[float]]:
    result: Dict[str, List[float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            raw = (row.get("dna_embedding") or row.get("dna_lm_embedding") or "").strip()
            if contig and raw and raw != ".":
                result[contig] = parse_vector(raw)
    dims = {len(value) for value in result.values()}
    if len(dims) > 1:
        raise ValueError("inconsistent DNA embedding dimensions")
    return result


def read_focus(path: Path) -> Set[str]:
    result: Set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader, None)
        if first and first[0].strip().lower() not in {"bin", "bin_id", "cluster"}:
            result.add(first[0].strip())
        for row in reader:
            if row and row[0].strip():
                result.add(row[0].strip())
    return result


def read_assignments(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        fields = list(reader.fieldnames)
        rows = [dict(row) for row in reader]
    mapping = {row["contig"]: row.get("bin", ".") for row in rows}
    return fields, rows, mapping


def prototype(z: torch.Tensor, indices: List[int]) -> torch.Tensor:
    return F.normalize(z[indices].mean(dim=0, keepdim=True), p=2, dim=-1)[0]


def train_local(
    matrix: torch.Tensor,
    left_idx: List[int],
    right_idx: List[int],
    args: argparse.Namespace,
) -> Tuple[ResidualMetricAdapter, Dict[str, float]]:
    dim = matrix.shape[1]
    adapter = ResidualMetricAdapter(dim, args.rank, args.residual_scale)
    optimizer = torch.optim.AdamW(adapter.parameters(), lr=args.learning_rate, weight_decay=1e-4)
    base = F.normalize(matrix, p=2, dim=-1)

    positive_pairs: List[Tuple[int, int]] = []
    for group in (left_idx, right_idx):
        for i, a in enumerate(group):
            for b in group[i + 1 :]:
                positive_pairs.append((a, b))
    negative_pairs = [(a, b) for a in left_idx for b in right_idx]
    if not negative_pairs:
        raise ValueError("local adaptation requires anchors on both sides")

    for _epoch in range(args.epochs):
        optimizer.zero_grad()
        z = adapter(matrix)
        if positive_pairs:
            pos = torch.stack([torch.dot(z[a], z[b]) for a, b in positive_pairs])
            pos_loss = (1.0 - pos).pow(2).mean()
        else:
            pos_loss = z.sum() * 0.0
        neg = torch.stack([torch.dot(z[a], z[b]) for a, b in negative_pairs])
        neg_loss = torch.relu(neg - args.negative_margin).pow(2).mean()
        drift = (1.0 - (z * base).sum(dim=-1)).pow(2).mean()
        loss = pos_loss + args.negative_weight * neg_loss + args.drift_weight * drift
        loss.backward()
        optimizer.step()

    adapter.eval()
    with torch.no_grad():
        z = adapter(matrix)
        left_proto = prototype(z, left_idx)
        right_proto = prototype(z, right_idx)
        cross = float(torch.dot(left_proto, right_proto).item())
        within_values = []
        for group in (left_idx, right_idx):
            for i, a in enumerate(group):
                for b in group[i + 1 :]:
                    within_values.append(float(torch.dot(z[a], z[b]).item()))
    return adapter, {
        "anchor_cross_prototype": cross,
        "anchor_within_mean": sum(within_values) / len(within_values) if within_values else float("nan"),
        "positive_anchor_pairs": len(positive_pairs),
        "negative_anchor_pairs": len(negative_pairs),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.rank < 1 or args.epochs < 1:
        raise SystemExit("--rank and --epochs must be positive")
    torch.manual_seed(args.seed)
    features = read_features(args.features)
    fields, rows, assignments = read_assignments(args.assignments)
    focus = read_focus(args.focus_bins)
    stage1 = json.loads(args.stage1_report.read_text(encoding="utf-8"))
    partitions = {
        entry["bin"]: entry
        for entry in stage1.get("bins", [])
        if entry.get("left_anchors") and entry.get("right_anchors")
    }
    members: Dict[str, List[str]] = defaultdict(list)
    for contig, bin_name in assignments.items():
        if bin_name not in {"", ".", "NA", "unbinned"} and contig in features:
            members[bin_name].append(contig)

    rewrites: Dict[str, str] = {}
    reports = []
    for bin_name in sorted(focus):
        entry = partitions.get(bin_name)
        current_members = sorted(members.get(bin_name, []))
        if entry is None or not current_members:
            continue
        left_anchors = [name for name in entry["left_anchors"] if name in features and name in current_members]
        right_anchors = [name for name in entry["right_anchors"] if name in features and name in current_members]
        if not left_anchors or not right_anchors:
            continue

        names = current_members
        index = {name: i for i, name in enumerate(names)}
        matrix = torch.tensor([features[name] for name in names], dtype=torch.float32)
        adapter, diagnostics = train_local(
            matrix,
            [index[name] for name in left_anchors],
            [index[name] for name in right_anchors],
            args,
        )
        with torch.no_grad():
            z = adapter(matrix)
            left_proto = prototype(z, [index[name] for name in left_anchors])
            right_proto = prototype(z, [index[name] for name in right_anchors])
        left_name = f"{bin_name}__bioA"
        right_name = f"{bin_name}__bioB"
        count_left = count_right = ambiguous = 0
        margins: List[float] = []
        left_anchor_set = set(left_anchors)
        right_anchor_set = set(right_anchors)
        for name in names:
            if name in left_anchor_set:
                rewrites[name] = left_name
                count_left += 1
                continue
            if name in right_anchor_set:
                rewrites[name] = right_name
                count_right += 1
                continue
            li = float(torch.dot(z[index[name]], left_proto).item())
            ri = float(torch.dot(z[index[name]], right_proto).item())
            margin = abs(li - ri)
            margins.append(margin)
            if margin < args.member_margin:
                # Abstain rather than inventing contamination. Keep it in the original bin.
                rewrites[name] = bin_name
                ambiguous += 1
            elif li >= ri:
                rewrites[name] = left_name
                count_left += 1
            else:
                rewrites[name] = right_name
                count_right += 1
        reports.append(
            {
                "bin": bin_name,
                "members": len(names),
                "left_anchors": len(left_anchors),
                "right_anchors": len(right_anchors),
                "assigned_left": count_left,
                "assigned_right": count_right,
                "ambiguous_kept": ambiguous,
                "member_margin_median": sorted(margins)[len(margins) // 2] if margins else 0.0,
                **diagnostics,
            }
        )

    for row in rows:
        contig = row.get("contig", "")
        if contig in rewrites:
            row["bin"] = rewrites[contig]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {"focused_bins": sorted(focus), "adapted_bins": len(reports), "member_margin": args.member_margin, "bins": reports}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"bridgebin-local-metric: focus={len(focus)} adapted={len(reports)} rewrites={len(rewrites)}")
    for item in reports:
        print(
            f"  {item['bin']} members={item['members']} anchors={item['left_anchors']}/{item['right_anchors']} "
            f"assigned={item['assigned_left']}/{item['assigned_right']} ambiguous={item['ambiguous_kept']} "
            f"cross={item['anchor_cross_prototype']:.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
