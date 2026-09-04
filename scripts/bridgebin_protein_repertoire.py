#!/usr/bin/env python3
"""Convert per-protein embeddings into genome-informative protein-repertoire vectors.

A naïve mean over all ESM-C proteins is dominated by conserved housekeeping proteins and
can make close genomes look *more* similar. BridgeBin instead learns an unsupervised
spherical-k-means codebook of protein embedding prototypes, assigns each ORF to a
prototype, and represents a contig by an L2-normalized TF-IDF histogram of prototypes.
Rare/accessory protein content therefore contributes more than ubiquitous proteins.

Input is the TSV produced by ``bridgebin_esmc_embed.py``:
    contig, protein_id, embedding, confidence, model

Commands:
  fit-transform  learn a codebook and emit contig repertoire features
  transform      apply an existing NPZ codebook to another sample

Requires numpy only; no sklearn dependency.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    fit = sub.add_parser("fit-transform")
    fit.add_argument("--embeddings", type=Path, required=True)
    fit.add_argument("--output", type=Path, required=True)
    fit.add_argument("--codebook-out", type=Path, required=True)
    fit.add_argument("--clusters", type=int, default=256)
    fit.add_argument("--max-train-proteins", type=int, default=20000)
    fit.add_argument("--iterations", type=int, default=20)
    fit.add_argument("--seed", type=int, default=43)

    transform = sub.add_parser("transform")
    transform.add_argument("--embeddings", type=Path, required=True)
    transform.add_argument("--output", type=Path, required=True)
    transform.add_argument("--codebook", type=Path, required=True)
    return parser.parse_args(argv)


def read_rows(path: Path) -> Iterator[Tuple[str, str, np.ndarray, float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        for row in reader:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            protein_id = (row.get("protein_id") or row.get("protein") or "").strip()
            raw = (row.get("embedding") or row.get("esm_embedding") or row.get("esmc_embedding") or "").strip()
            if not contig or not raw:
                continue
            vector = np.fromstring(raw, sep=",", dtype=np.float32)
            if vector.size == 0 or not np.isfinite(vector).all():
                continue
            confidence_raw = (row.get("confidence") or row.get("score") or "1").strip()
            confidence = float(confidence_raw) if confidence_raw not in {"", ".", "NA"} else 1.0
            if not math.isfinite(confidence):
                continue
            yield contig, protein_id, vector, max(0.0, min(1.0, confidence))


def normalize_rows(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    return matrix / norms


def reservoir_sample(path: Path, limit: int, seed: int) -> np.ndarray:
    if limit < 1:
        raise ValueError("--max-train-proteins must be positive")
    rng = random.Random(seed)
    sample: List[np.ndarray] = []
    dimension: Optional[int] = None
    seen = 0
    for _contig, _protein, vector, confidence in read_rows(path):
        if confidence <= 0.0:
            continue
        if dimension is None:
            dimension = int(vector.size)
        if int(vector.size) != dimension:
            raise ValueError(f"inconsistent protein embedding dimension: {vector.size} vs {dimension}")
        seen += 1
        if len(sample) < limit:
            sample.append(vector.copy())
        else:
            replacement = rng.randrange(seen)
            if replacement < limit:
                sample[replacement] = vector.copy()
    if not sample:
        raise ValueError("no usable protein embeddings")
    matrix = np.stack(sample).astype(np.float32, copy=False)
    return normalize_rows(matrix)


def kmeans_plus_plus(data: np.ndarray, clusters: int, seed: int) -> np.ndarray:
    n, dimension = data.shape
    if clusters < 2 or clusters > n:
        raise ValueError(f"--clusters must be in [2,{n}]")
    rng = np.random.default_rng(seed)
    centers = np.empty((clusters, dimension), dtype=np.float32)
    first = int(rng.integers(n))
    centers[0] = data[first]
    best_similarity = data @ centers[0]
    for index in range(1, clusters):
        distance = np.maximum(0.0, 1.0 - best_similarity)
        weights = distance * distance
        total = float(weights.sum())
        if total <= 1e-12:
            chosen = int(rng.integers(n))
        else:
            chosen = int(rng.choice(n, p=weights / total))
        centers[index] = data[chosen]
        best_similarity = np.maximum(best_similarity, data @ centers[index])
    return normalize_rows(centers)


def spherical_kmeans(
    data: np.ndarray, clusters: int, iterations: int, seed: int
) -> Tuple[np.ndarray, Dict[str, float]]:
    if iterations < 1:
        raise ValueError("--iterations must be positive")
    centers = kmeans_plus_plus(data, clusters, seed)
    rng = np.random.default_rng(seed + 1)
    previous = None
    mean_similarity = 0.0
    for iteration in range(iterations):
        similarities = data @ centers.T
        labels = np.argmax(similarities, axis=1)
        best = similarities[np.arange(data.shape[0]), labels]
        mean_similarity = float(best.mean())
        if previous is not None and np.array_equal(labels, previous):
            break
        previous = labels.copy()
        updated = np.zeros_like(centers)
        counts = np.bincount(labels, minlength=clusters)
        np.add.at(updated, labels, data)
        for cluster in range(clusters):
            if counts[cluster] == 0:
                updated[cluster] = data[int(rng.integers(data.shape[0]))]
        centers = normalize_rows(updated)
    return centers, {
        "train_proteins": float(data.shape[0]),
        "dimension": float(data.shape[1]),
        "clusters": float(clusters),
        "mean_train_cosine": mean_similarity,
    }


def save_codebook(path: Path, centers: np.ndarray, metadata: Dict[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, centers=centers.astype(np.float32), metadata=json.dumps(metadata))


def load_codebook(path: Path) -> np.ndarray:
    with np.load(path, allow_pickle=False) as archive:
        centers = np.asarray(archive["centers"], dtype=np.float32)
    if centers.ndim != 2 or centers.shape[0] < 2:
        raise ValueError("invalid protein repertoire codebook")
    return normalize_rows(centers)


def assign_repertoire(
    path: Path, centers: np.ndarray
) -> Tuple[Dict[str, np.ndarray], Dict[str, Tuple[float, int]], int]:
    clusters, dimension = centers.shape
    histograms: Dict[str, np.ndarray] = {}
    quality: Dict[str, Tuple[float, int]] = defaultdict(lambda: (0.0, 0))
    total = 0
    for contig, _protein_id, vector, confidence in read_rows(path):
        if vector.size != dimension:
            raise ValueError(
                f"embedding dimension {vector.size} does not match codebook dimension {dimension}"
            )
        norm = float(np.linalg.norm(vector))
        if norm <= 1e-12:
            continue
        normalized = vector / norm
        similarities = centers @ normalized
        cluster = int(np.argmax(similarities))
        histogram = histograms.setdefault(contig, np.zeros(clusters, dtype=np.float32))
        # Downweight proteins that are both low confidence and poorly represented by
        # the codebook; negative cosine cannot add evidence.
        match = max(0.0, float(similarities[cluster]))
        histogram[cluster] += confidence * (0.5 + 0.5 * match)
        quality_sum, count = quality[contig]
        quality[contig] = (quality_sum + confidence * match, count + 1)
        total += 1
    return histograms, quality, total


def tfidf(histograms: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    if not histograms:
        return {}
    matrix = np.stack(list(histograms.values()))
    document_frequency = (matrix > 0).sum(axis=0)
    n_documents = matrix.shape[0]
    idf = np.log((1.0 + n_documents) / (1.0 + document_frequency)) + 1.0
    result: Dict[str, np.ndarray] = {}
    for contig, histogram in histograms.items():
        weighted = np.log1p(histogram) * idf
        norm = float(np.linalg.norm(weighted))
        if norm > 1e-12:
            weighted = weighted / norm
        result[contig] = weighted.astype(np.float32, copy=False)
    return result


def write_features(
    output: Path,
    histograms: Dict[str, np.ndarray],
    quality: Dict[str, Tuple[float, int]],
    total_proteins: int,
) -> None:
    features = tfidf(histograms)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "protein_repertoire", "repertoire_confidence", "protein_count"])
        for contig in sorted(features):
            quality_sum, count = quality[contig]
            mean_quality = quality_sum / max(1, count)
            amount = 1.0 - math.exp(-count / 8.0)
            confidence = max(0.0, min(1.0, mean_quality * amount))
            writer.writerow(
                [
                    contig,
                    ",".join(f"{float(value):.7g}" for value in features[contig]),
                    f"{confidence:.6f}",
                    count,
                ]
            )
    print(
        f"bridgebin-protein-repertoire: contigs={len(features)} proteins={total_proteins} output={output}"
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.command == "fit-transform":
        data = reservoir_sample(args.embeddings, args.max_train_proteins, args.seed)
        centers, metadata = spherical_kmeans(data, args.clusters, args.iterations, args.seed)
        save_codebook(args.codebook_out, centers, metadata)
        histograms, quality, total = assign_repertoire(args.embeddings, centers)
        write_features(args.output, histograms, quality, total)
        print("bridgebin-protein-repertoire-codebook:", json.dumps(metadata, sort_keys=True))
        return 0
    if args.command == "transform":
        centers = load_codebook(args.codebook)
        histograms, quality, total = assign_repertoire(args.embeddings, centers)
        write_features(args.output, histograms, quality, total)
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
