#!/usr/bin/env python3
"""Pool DNA-LM windows for genome identity without mixing dispersion into cosine space.

BridgeBin's generic feature builder keeps mean and standard deviation together because that
is useful for heterogeneous multimodal features.  For DNABERT-S genome identity, however,
the foundation-model embedding itself is the identity coordinate: concatenating per-window
standard deviation changes the geometry, and multiplying cosine similarity by window
confidence can change pair ranking.

This small adapter therefore emits the same BridgeBin feature TSV schema while using only
the confidence-weighted mean DNA embedding.  Window confidence is retained as metadata for
future uncertainty handling, but it is not concatenated into the identity vector.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


@dataclass
class Pool:
    weighted_sum: List[float] = field(default_factory=list)
    total_weight: float = 0.0
    observations: int = 0

    def add(self, vector: Sequence[float], weight: float) -> None:
        if not vector or weight <= 0.0:
            return
        if not self.weighted_sum:
            self.weighted_sum = [0.0] * len(vector)
        if len(vector) != len(self.weighted_sum):
            raise ValueError(
                f"inconsistent DNA embedding dimension: expected {len(self.weighted_sum)}, got {len(vector)}"
            )
        for index, value in enumerate(vector):
            self.weighted_sum[index] += weight * value
        self.total_weight += weight
        self.observations += 1

    def mean(self) -> List[float]:
        if self.total_weight <= 0.0:
            return []
        return [value / self.total_weight for value in self.weighted_sum]

    def confidence(self, saturation: float = 2.0) -> float:
        if self.total_weight <= 0.0 or self.observations <= 0:
            return 0.0
        mean_weight = min(1.0, self.total_weight / self.observations)
        amount = 1.0 - math.exp(-self.observations / saturation)
        return max(0.0, min(1.0, mean_weight * amount))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dna-embeddings", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


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
    values = [float(value) for value in raw.split(",") if value.strip()]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("DNA embedding must contain finite comma-separated floats")
    return values


def parse_probability(raw: str, default: float = 1.0) -> float:
    if not raw:
        return default
    value = float(raw)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"confidence outside [0,1]: {raw!r}")
    return value


def format_vector(values: Sequence[float]) -> str:
    return ",".join(f"{value:.7g}" for value in values)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    pooled: Dict[str, Pool] = defaultdict(Pool)
    windows = 0
    for row in rows(args.dna_embeddings):
        contig = first(row, ("contig", "contig_id", "sequence"))
        raw = first(row, ("embedding", "dna_embedding", "dna_lm_embedding", "dnabert_embedding"))
        if not contig or not raw:
            continue
        weight = parse_probability(first(row, ("confidence", "score")), 1.0)
        pooled[contig].add(parse_vector(raw), weight)
        windows += 1

    if not pooled:
        raise SystemExit("no usable DNA embeddings")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "contig",
                "taxonomy",
                "taxonomy_confidence",
                "dna_embedding",
                "dna_confidence",
                "gene_profile",
                "gene_confidence",
                "esm_embedding",
                "protein_confidence",
            ]
        )
        for contig in sorted(pooled):
            feature = pooled[contig]
            writer.writerow(
                [
                    contig,
                    ".",
                    "0",
                    format_vector(feature.mean()),
                    f"{feature.confidence():.6f}",
                    ".",
                    "0",
                    ".",
                    "0",
                ]
            )

    print(
        f"bridgebin-dna-identity-features: contigs={len(pooled)} windows={windows} pooling=weighted-mean"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
