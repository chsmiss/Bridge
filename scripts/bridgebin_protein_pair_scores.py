#!/usr/bin/env python3
"""Compute contig-pair protein-set evidence from per-protein language-model embeddings.

A single mean protein embedding is a poor genome barcode: one conserved housekeeping
protein can dominate it.  This scorer keeps the ORF set structure.  For every candidate
contig pair it computes bidirectional best matches between proteins and summarizes the
*distribution* of those matches.  Strong evidence therefore requires several proteins on
both contigs to agree, not just one universal homolog.

The output preserves the input candidate table and adds ``protein_set_similarity``,
``protein_set_coverage`` and ``protein_confidence``.  These are soft biological evidence;
absence of a match must not be treated as a hard cannot-link for short contigs.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--protein-embeddings", type=Path, required=True)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--match-threshold", type=float, default=0.55)
    return p.parse_args(argv)


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        yield from reader


def vector(raw: str) -> List[float]:
    values = [float(piece) for piece in raw.split(",") if piece.strip()]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("protein embedding must contain finite comma-separated floats")
    norm = math.sqrt(sum(value * value for value in values))
    if norm <= 1e-12:
        raise ValueError("zero-norm protein embedding")
    return [value / norm for value in values]


def cosine(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        return 0.0
    return max(0.0, min(1.0, sum(a * b for a, b in zip(left, right))))


def median(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def read_embeddings(path: Path) -> Dict[str, List[Tuple[List[float], float]]]:
    result: Dict[str, List[Tuple[List[float], float]]] = defaultdict(list)
    expected_dim: Optional[int] = None
    for row in rows(path):
        contig = (row.get("contig") or row.get("contig_id") or "").strip()
        raw = (row.get("embedding") or row.get("esm_embedding") or row.get("esmc_embedding") or "").strip()
        if not contig or not raw or raw == ".":
            continue
        embedding = vector(raw)
        if expected_dim is None:
            expected_dim = len(embedding)
        elif len(embedding) != expected_dim:
            raise ValueError("inconsistent protein embedding dimensions")
        confidence = float((row.get("confidence") or "1").strip() or 1.0)
        confidence = max(0.0, min(1.0, confidence))
        result[contig].append((embedding, confidence))
    return result


def pair_score(
    left: Sequence[Tuple[List[float], float]],
    right: Sequence[Tuple[List[float], float]],
    threshold: float,
) -> Tuple[Optional[float], float, float]:
    if not left or not right:
        return None, 0.0, 0.0
    matrix = [[cosine(a, b) for b, _ in right] for a, _ in left]
    left_best = [max(row) for row in matrix]
    right_best = [max(matrix[i][j] for i in range(len(matrix))) for j in range(len(right))]
    robust_similarity = 0.5 * (median(left_best) + median(right_best))
    left_cov = sum(value >= threshold for value in left_best) / len(left_best)
    right_cov = sum(value >= threshold for value in right_best) / len(right_best)
    coverage = math.sqrt(left_cov * right_cov)
    # Require both similarity and breadth.  A lone conserved protein should be weak evidence.
    score = robust_similarity * (0.35 + 0.65 * coverage)
    left_conf = sum(conf for _vec, conf in left) / len(left)
    right_conf = sum(conf for _vec, conf in right) / len(right)
    amount = 1.0 - math.exp(-min(len(left), len(right)) / 2.0)
    confidence = min(left_conf, right_conf) * amount
    return score, coverage, max(0.0, min(1.0, confidence))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.match_threshold <= 1.0:
        raise SystemExit("--match-threshold must be in [0,1]")
    embeddings = read_embeddings(args.protein_embeddings)
    with args.pairs.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{args.pairs}: missing header")
        fields = list(reader.fieldnames)
        for extra in ("protein_set_similarity", "protein_set_coverage", "protein_confidence"):
            if extra not in fields:
                fields.append(extra)
        pair_rows = [dict(row) for row in reader]

    scored = 0
    for row in pair_rows:
        left = (row.get("left") or row.get("source") or "").strip()
        right = (row.get("right") or row.get("target") or "").strip()
        score, coverage, confidence = pair_score(
            embeddings.get(left, []), embeddings.get(right, []), args.match_threshold
        )
        row["protein_set_similarity"] = "." if score is None else f"{score:.8f}"
        row["protein_set_coverage"] = f"{coverage:.8f}"
        row["protein_confidence"] = f"{confidence:.8f}"
        if score is not None:
            scored += 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pair_rows)
    print(f"bridgebin-protein-pairs: candidates={len(pair_rows)} scored={scored} contigs={len(embeddings)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
