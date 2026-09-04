#!/usr/bin/env python3
"""Parameter-efficient metric transfer for BridgeBin DNA foundation-model embeddings.

The DNABERT-S backbone stays frozen.  We train a small residual adapter on top of
contig-level DNA embeddings so that unseen genomes remain separated in the embedding
space without teaching the model benchmark-specific absolute probabilities.

The training objective is deliberately asymmetric for binning:

* same-genome pairs are pulled together;
* different-genome pairs are pushed below a cosine margin;
* negatives that were already very similar in the frozen representation receive larger
  weight (hard-negative curriculum);
* a small drift penalty keeps the adapted representation near the pretrained geometry.

Validation holds out complete genome groups when group columns are present.  ``transform``
preserves the BridgeBin feature-table schema and only replaces ``dna_embedding``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


class ResidualMetricAdapter(nn.Module):
    def __init__(self, dim: int, rank: int, residual_scale: float) -> None:
        super().__init__()
        self.dim = dim
        self.rank = rank
        self.residual_scale = residual_scale
        self.down = nn.Linear(dim, rank, bias=False)
        self.up = nn.Linear(rank, dim, bias=False)
        nn.init.normal_(self.down.weight, mean=0.0, std=0.02)
        nn.init.zeros_(self.up.weight)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        base = F.normalize(x, p=2, dim=-1)
        delta = self.up(F.gelu(self.down(base)))
        return F.normalize(base + self.residual_scale * delta, p=2, dim=-1)


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        yield from reader


def first(row: Dict[str, str], names: Sequence[str], default: str = "") -> str:
    for name in names:
        value = row.get(name)
        if value is not None and value.strip() and value.strip() not in {".", "NA", "na"}:
            return value.strip()
    return default


def parse_vector(raw: str) -> List[float]:
    values = [float(piece) for piece in raw.split(",") if piece.strip()]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("DNA embedding must contain finite comma-separated floats")
    return values


def read_feature_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]], Dict[str, List[float]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "contig" not in reader.fieldnames:
            raise ValueError(f"{path}: feature table needs a contig column")
        fields = list(reader.fieldnames)
        output_rows = [dict(row) for row in reader]
    vectors: Dict[str, List[float]] = {}
    for row in output_rows:
        contig = first(row, ("contig", "contig_id"))
        raw = first(row, ("dna_embedding", "dna_lm_embedding", "dnabert_embedding"))
        if contig and raw:
            vectors[contig] = parse_vector(raw)
    if not vectors:
        raise ValueError(f"{path}: no DNA embeddings")
    dimensions = {len(vector) for vector in vectors.values()}
    if len(dimensions) != 1:
        raise ValueError(f"{path}: inconsistent DNA embedding dimensions: {sorted(dimensions)}")
    return fields, output_rows, vectors


def hash_unit(value: str, seed: int) -> float:
    digest = hashlib.blake2b(f"{value}\0{seed}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") / float(2**64 - 1)


def read_pairs(path: Path, vectors: Dict[str, List[float]]) -> List[Tuple[str, str, int, str, str]]:
    pairs: List[Tuple[str, str, int, str, str]] = []
    for row in rows(path):
        left = first(row, ("left", "source", "contig_a", "contig1"))
        right = first(row, ("right", "target", "contig_b", "contig2"))
        label_raw = first(row, ("label", "same_genome", "truth"))
        if not left or not right or left == right or left not in vectors or right not in vectors:
            continue
        if not label_raw:
            continue
        label_float = float(label_raw)
        if label_float not in {0.0, 1.0}:
            raise ValueError(f"pair label must be 0/1, got {label_raw!r}")
        left_group = first(row, ("left_group", "left_genome", "genome_left", "source_genome"))
        right_group = first(row, ("right_group", "right_genome", "genome_right", "target_genome"))
        pairs.append((left, right, int(label_float), left_group, right_group))
    if not pairs:
        raise ValueError("no usable labelled pairs")
    return pairs


def split_pairs(
    pairs: Sequence[Tuple[str, str, int, str, str]],
    validation_fraction: float,
    seed: int,
    require_group_holdout: bool,
) -> Tuple[List[int], List[int], Dict[str, object]]:
    groups = sorted({group for *_prefix, lg, rg in pairs for group in (lg, rg) if group})
    complete_groups = bool(groups) and all(lg and rg for *_prefix, lg, rg in pairs)
    if complete_groups and len(groups) >= 3:
        ordered = sorted(groups, key=lambda group: (hash_unit(group, seed), group))
        holdout_count = max(1, min(len(ordered) - 2, round(len(ordered) * validation_fraction)))
        held = set(ordered[:holdout_count])
        train_idx = [i for i, (*_head, lg, rg) in enumerate(pairs) if lg not in held and rg not in held]
        valid_idx = [i for i, (*_head, lg, rg) in enumerate(pairs) if lg in held or rg in held]
        protocol = {
            "split_protocol": "whole_genome_holdout",
            "validation_groups": sorted(held),
            "all_groups": len(groups),
        }
    else:
        if require_group_holdout:
            raise ValueError("--require-group-holdout needs complete group metadata and >=3 groups")
        train_idx = []
        valid_idx = []
        for i, (left, right, *_rest) in enumerate(pairs):
            key = "\0".join(sorted((left, right)))
            (valid_idx if hash_unit(key, seed) < validation_fraction else train_idx).append(i)
        protocol = {"split_protocol": "pair_hash_fallback", "validation_groups": [], "all_groups": len(groups)}
    if not train_idx or not valid_idx:
        raise ValueError("empty train/validation split")
    return train_idx, valid_idx, protocol


def roc_auc(scores: Sequence[Tuple[float, int]]) -> float:
    positives = [score for score, label in scores if label == 1]
    negatives = [score for score, label in scores if label == 0]
    if not positives or not negatives:
        return float("nan")
    wins = 0.0
    for pos in positives:
        for neg in negatives:
            if pos > neg:
                wins += 1.0
            elif pos == neg:
                wins += 0.5
    return wins / (len(positives) * len(negatives))


def percentile(values: Sequence[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * q
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return ordered[lo]
    weight = position - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def evaluate(
    adapter: ResidualMetricAdapter,
    matrix: torch.Tensor,
    pair_index: Dict[str, int],
    pairs: Sequence[Tuple[str, str, int, str, str]],
    indices: Sequence[int],
) -> Dict[str, float]:
    adapter.eval()
    with torch.no_grad():
        z = adapter(matrix)
    scored: List[Tuple[float, int]] = []
    for pair_i in indices:
        left, right, label, _lg, _rg = pairs[pair_i]
        score = float(torch.dot(z[pair_index[left]], z[pair_index[right]]).item())
        scored.append((score, label))
    pos = [score for score, label in scored if label == 1]
    neg = [score for score, label in scored if label == 0]
    return {
        "auc": roc_auc(scored),
        "positive_median": percentile(pos, 0.5),
        "negative_median": percentile(neg, 0.5),
        "positive_q10": percentile(pos, 0.10),
        "negative_q90": percentile(neg, 0.90),
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    train = sub.add_parser("train", help="train residual metric adapter on labelled genome pairs")
    train.add_argument("--features", type=Path, required=True)
    train.add_argument("--pairs", type=Path, required=True)
    train.add_argument("--model-out", type=Path, required=True)
    train.add_argument("--validation-fraction", type=float, default=0.28)
    train.add_argument("--rank", type=int, default=64)
    train.add_argument("--residual-scale", type=float, default=0.35)
    train.add_argument("--epochs", type=int, default=500)
    train.add_argument("--learning-rate", type=float, default=3e-3)
    train.add_argument("--weight-decay", type=float, default=1e-4)
    train.add_argument("--negative-margin", type=float, default=0.25)
    train.add_argument("--negative-weight", type=float, default=3.0)
    train.add_argument("--hard-negative-scale", type=float, default=5.0)
    train.add_argument("--drift-weight", type=float, default=0.10)
    train.add_argument("--seed", type=int, default=43)
    train.add_argument("--require-group-holdout", action="store_true")

    transform = sub.add_parser("transform", help="apply a trained adapter to a BridgeBin feature TSV")
    transform.add_argument("--features", type=Path, required=True)
    transform.add_argument("--model", type=Path, required=True)
    transform.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def train(args: argparse.Namespace) -> int:
    if not 0.0 < args.validation_fraction < 0.5:
        raise ValueError("--validation-fraction must be in (0, 0.5)")
    if args.rank < 1 or args.epochs < 1:
        raise ValueError("--rank and --epochs must be positive")
    if args.negative_margin < -1.0 or args.negative_margin > 1.0:
        raise ValueError("--negative-margin must be in [-1,1]")

    torch.manual_seed(args.seed)
    random.seed(args.seed)
    _fields, _rows, vectors = read_feature_rows(args.features)
    pairs = read_pairs(args.pairs, vectors)
    train_idx, valid_idx, protocol = split_pairs(
        pairs, args.validation_fraction, args.seed, args.require_group_holdout
    )
    if {pairs[i][2] for i in train_idx} != {0, 1}:
        raise ValueError("training split must contain positive and negative pairs")
    if {pairs[i][2] for i in valid_idx} != {0, 1}:
        raise ValueError("validation split must contain positive and negative pairs")

    names = sorted(vectors)
    index = {name: i for i, name in enumerate(names)}
    matrix = torch.tensor([vectors[name] for name in names], dtype=torch.float32)
    dim = matrix.shape[1]
    adapter = ResidualMetricAdapter(dim, args.rank, args.residual_scale)
    optimizer = torch.optim.AdamW(
        adapter.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay
    )

    base = F.normalize(matrix, p=2, dim=-1)
    pair_left = torch.tensor([index[pairs[i][0]] for i in train_idx], dtype=torch.long)
    pair_right = torch.tensor([index[pairs[i][1]] for i in train_idx], dtype=torch.long)
    labels = torch.tensor([pairs[i][2] for i in train_idx], dtype=torch.float32)
    base_cos = (base[pair_left] * base[pair_right]).sum(dim=-1).detach()

    best_state = None
    best_auc = -1.0
    best_epoch = -1
    for epoch in range(args.epochs):
        adapter.train()
        optimizer.zero_grad()
        z = adapter(matrix)
        cos = (z[pair_left] * z[pair_right]).sum(dim=-1)
        pos_mask = labels > 0.5
        neg_mask = ~pos_mask
        pos_loss = (1.0 - cos[pos_mask]).pow(2).mean() if pos_mask.any() else cos.sum() * 0.0
        if neg_mask.any():
            hardness = 1.0 + args.hard_negative_scale * torch.relu(base_cos[neg_mask] - args.negative_margin)
            neg_term = torch.relu(cos[neg_mask] - args.negative_margin).pow(2)
            neg_loss = (hardness * neg_term).mean()
        else:
            neg_loss = cos.sum() * 0.0
        drift = (1.0 - (z * base).sum(dim=-1)).pow(2).mean()
        loss = pos_loss + args.negative_weight * neg_loss + args.drift_weight * drift
        loss.backward()
        optimizer.step()

        if epoch % 10 == 0 or epoch + 1 == args.epochs:
            metrics = evaluate(adapter, matrix, index, pairs, valid_idx)
            auc = metrics["auc"]
            if math.isfinite(auc) and auc > best_auc + 1e-9:
                best_auc = auc
                best_epoch = epoch
                best_state = {key: value.detach().cpu().clone() for key, value in adapter.state_dict().items()}

    if best_state is not None:
        adapter.load_state_dict(best_state)
    train_metrics = evaluate(adapter, matrix, index, pairs, train_idx)
    valid_metrics = evaluate(adapter, matrix, index, pairs, valid_idx)
    base_adapter = ResidualMetricAdapter(dim, args.rank, 0.0)
    base_metrics = evaluate(base_adapter, matrix, index, pairs, valid_idx)

    payload = {
        "version": 1,
        "input_dim": dim,
        "rank": args.rank,
        "residual_scale": args.residual_scale,
        "state_dict": adapter.state_dict(),
        "metadata": {
            "best_epoch": best_epoch,
            "train_pairs": len(train_idx),
            "validation_pairs": len(valid_idx),
            "negative_margin": args.negative_margin,
            "negative_weight": args.negative_weight,
            "hard_negative_scale": args.hard_negative_scale,
            "drift_weight": args.drift_weight,
            "base_validation": base_metrics,
            "adapted_train": train_metrics,
            "adapted_validation": valid_metrics,
            **protocol,
        },
    }
    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(payload, args.model_out)
    print(json.dumps(payload["metadata"], sort_keys=True))
    return 0


def transform(args: argparse.Namespace) -> int:
    fields, table_rows, vectors = read_feature_rows(args.features)
    payload = torch.load(args.model, map_location="cpu", weights_only=False)
    dim = int(payload["input_dim"])
    rank = int(payload["rank"])
    residual_scale = float(payload["residual_scale"])
    adapter = ResidualMetricAdapter(dim, rank, residual_scale)
    adapter.load_state_dict(payload["state_dict"])
    adapter.eval()

    names = sorted(vectors)
    if any(len(vectors[name]) != dim for name in names):
        raise ValueError("feature DNA dimension does not match metric-transfer model")
    matrix = torch.tensor([vectors[name] for name in names], dtype=torch.float32)
    with torch.no_grad():
        adapted = adapter(matrix).cpu().tolist()
    mapped = {
        name: ",".join(f"{value:.7g}" for value in vector)
        for name, vector in zip(names, adapted)
    }
    for row in table_rows:
        contig = first(row, ("contig", "contig_id"))
        if contig in mapped:
            row["dna_embedding"] = mapped[contig]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table_rows)
    print(f"bridgebin-metric-transfer: transformed={len(mapped)} output={args.output}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.command == "train":
        return train(args)
    if args.command == "transform":
        return transform(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
