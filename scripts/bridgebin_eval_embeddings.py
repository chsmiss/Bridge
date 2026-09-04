#!/usr/bin/env python3
"""Evaluate contig embedding separability against benchmark truth.

This is a benchmark-only diagnostic; truth is never consumed by the production pipeline.
It quantifies whether a pretrained representation actually supplies genome identity beyond
coverage/composition before we spend effort integrating it into the binner.

Reported metrics include pairwise ROC AUC, nearest-neighbor genome accuracy, same/different
cosine quantiles, and precision/recall at a requested high-precision threshold. Optional
``--focus A,B`` reports a specific hard-negative genome pair such as E. coli/Salmonella.
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


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--column", default="dna_embedding")
    parser.add_argument("--max-contigs-per-genome", type=int, default=80)
    parser.add_argument("--precision-target", type=float, default=0.99)
    parser.add_argument("--focus", action="append", default=[], help="genomeA,genomeB")
    parser.add_argument("--seed", type=int, default=43)
    return parser.parse_args(argv)


def parse_vector(raw: str) -> List[float]:
    values = [float(piece) for piece in raw.split(",") if piece.strip()]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("invalid embedding vector")
    norm = math.sqrt(sum(value * value for value in values))
    if norm <= 1e-12:
        raise ValueError("zero-norm embedding vector")
    return [value / norm for value in values]


def read_features(path: Path, column: str) -> Dict[str, List[float]]:
    result: Dict[str, List[float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing embedding column {column!r}")
        id_column = next(
            (name for name in ("contig", "contig_id", "sequence") if name in reader.fieldnames),
            None,
        )
        if id_column is None:
            raise ValueError(f"{path}: missing contig identifier column")
        dimension: Optional[int] = None
        for row in reader:
            contig = (row.get(id_column) or "").strip()
            raw = (row.get(column) or "").strip()
            if not contig or raw in {"", ".", "NA", "na"}:
                continue
            vector = parse_vector(raw)
            if dimension is None:
                dimension = len(vector)
            if len(vector) != dimension:
                raise ValueError(f"inconsistent embedding dimension for {contig}")
            result[contig] = vector
    return result


def read_truth(path: Path) -> Dict[str, Tuple[str, int]]:
    result: Dict[str, Tuple[str, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            genome = (row.get("genome") or row.get("species") or "").strip()
            if not contig or not genome:
                continue
            length = int(float(row.get("length") or 0))
            eligible = (row.get("eligible") or "1").strip().lower() not in {
                "0",
                "false",
                "no",
            }
            if eligible:
                result[contig] = (genome, length)
    return result


def dot(left: Sequence[float], right: Sequence[float]) -> float:
    return max(-1.0, min(1.0, sum(a * b for a, b in zip(left, right))))


def quantile(values: Sequence[float], fraction: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return ordered[lo]
    mix = position - lo
    return ordered[lo] * (1.0 - mix) + ordered[hi] * mix


def summarize(values: Sequence[float]) -> Dict[str, float]:
    return {
        "n": len(values),
        "q05": quantile(values, 0.05),
        "q25": quantile(values, 0.25),
        "median": quantile(values, 0.50),
        "q75": quantile(values, 0.75),
        "q95": quantile(values, 0.95),
        "mean": sum(values) / max(1, len(values)),
    }


def auc_mann_whitney(positive: Sequence[float], negative: Sequence[float]) -> float:
    if not positive or not negative:
        return float("nan")
    combined = [(score, 1) for score in positive] + [(score, 0) for score in negative]
    combined.sort(key=lambda item: item[0])
    rank_sum_positive = 0.0
    index = 0
    while index < len(combined):
        end = index + 1
        while end < len(combined) and combined[end][0] == combined[index][0]:
            end += 1
        average_rank = (index + 1 + end) / 2.0
        rank_sum_positive += average_rank * sum(label for _, label in combined[index:end])
        index = end
    n_pos = len(positive)
    n_neg = len(negative)
    u = rank_sum_positive - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def threshold_at_precision(
    positive: Sequence[float], negative: Sequence[float], target: float
) -> Tuple[float, float, float]:
    labelled = [(score, 1) for score in positive] + [(score, 0) for score in negative]
    labelled.sort(key=lambda item: item[0], reverse=True)
    tp = fp = 0
    total_positive = len(positive)
    chosen = 1.0
    chosen_precision = 1.0
    chosen_recall = 0.0
    index = 0
    while index < len(labelled):
        score = labelled[index][0]
        end = index
        while end < len(labelled) and labelled[end][0] == score:
            if labelled[end][1]:
                tp += 1
            else:
                fp += 1
            end += 1
        precision = tp / max(1, tp + fp)
        recall = tp / max(1, total_positive)
        if precision >= target and recall >= chosen_recall:
            chosen = score
            chosen_precision = precision
            chosen_recall = recall
        index = end
    return chosen, chosen_precision, chosen_recall


def select_contigs(
    features: Dict[str, List[float]],
    truth: Dict[str, Tuple[str, int]],
    max_per_genome: int,
    seed: int,
) -> List[str]:
    grouped: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    for contig in set(features).intersection(truth):
        genome, length = truth[contig]
        grouped[genome].append((contig, length))
    selected: List[str] = []
    for genome, values in sorted(grouped.items()):
        values.sort(key=lambda item: (-item[1], item[0]))
        if max_per_genome > 0 and len(values) > max_per_genome:
            # Keep half long representatives and half deterministic diversity sample.
            keep_long = max(1, max_per_genome // 2)
            head = values[:keep_long]
            tail = values[keep_long:]
            rng = random.Random(f"{seed}:{genome}")
            rng.shuffle(tail)
            values = head + tail[: max_per_genome - len(head)]
        selected.extend(contig for contig, _ in values)
    return selected


def evaluate(args: argparse.Namespace) -> Dict[str, object]:
    if not 0.5 < args.precision_target < 1.0:
        raise ValueError("--precision-target must be in (0.5,1)")
    features = read_features(args.features, args.column)
    truth = read_truth(args.truth)
    selected = select_contigs(
        features, truth, args.max_contigs_per_genome, args.seed
    )
    if len(selected) < 3:
        raise ValueError("too few benchmark contigs with embeddings and truth")

    same: List[float] = []
    different: List[float] = []
    focus_requested: List[Tuple[str, str]] = []
    for raw in args.focus:
        parts = [part.strip() for part in raw.split(",")]
        if len(parts) != 2 or not all(parts):
            raise ValueError(f"--focus must be genomeA,genomeB, got {raw!r}")
        focus_requested.append((parts[0], parts[1]))
    focus_values: Dict[str, List[float]] = {
        " vs ".join(pair): [] for pair in focus_requested
    }

    nearest_correct = 0
    for left_pos, left in enumerate(selected):
        left_genome = truth[left][0]
        best_score = -2.0
        best_genome = None
        for right_pos, right in enumerate(selected):
            if left_pos == right_pos:
                continue
            score = dot(features[left], features[right])
            right_genome = truth[right][0]
            if right_pos > left_pos:
                if left_genome == right_genome:
                    same.append(score)
                else:
                    different.append(score)
                for pair in focus_requested:
                    if {left_genome, right_genome} == set(pair):
                        focus_values[" vs ".join(pair)].append(score)
            if score > best_score:
                best_score = score
                best_genome = right_genome
        nearest_correct += int(best_genome == left_genome)

    threshold, precision, recall = threshold_at_precision(
        same, different, args.precision_target
    )
    genomes = sorted({truth[contig][0] for contig in selected})
    return {
        "column": args.column,
        "selected_contigs": len(selected),
        "genomes": genomes,
        "genome_count": len(genomes),
        "pair_auc": auc_mann_whitney(same, different),
        "nearest_neighbor_genome_accuracy": nearest_correct / len(selected),
        "same_genome_cosine": summarize(same),
        "different_genome_cosine": summarize(different),
        "high_precision_threshold": {
            "target": args.precision_target,
            "threshold": threshold,
            "precision": precision,
            "recall": recall,
        },
        "focus_cross_genome_cosine": {
            name: summarize(values) for name, values in focus_values.items()
        },
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    metrics = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(metrics, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
