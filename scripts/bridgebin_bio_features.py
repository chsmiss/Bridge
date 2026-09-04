#!/usr/bin/env python3
"""Aggregate DNA, gene annotation, and protein-language-model outputs for BridgeBin.

This utility deliberately separates expensive Python/GPU inference from the Rust
partitioner. It accepts generic contig DNA embeddings, per-gene annotations, per-protein
embeddings, and taxonomy, then produces a header-based feature TSV.

Gene profiles use deterministic feature hashing, so outputs from GENERanno, Prodigal +
HMM/DIAMOND annotation, or another gene annotator can share the same interface. Protein
embeddings may come from ESM-C or another protein model. DNA embeddings may come from
DNABERT-S, GENERanno-base, or another nucleotide foundation model. For variable numbers
of windows/ORFs we preserve more information than a naive mean by concatenating the
confidence-weighted mean and standard deviation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


@dataclass
class WeightedVectors:
    weighted_sum: List[float] = field(default_factory=list)
    weighted_sq_sum: List[float] = field(default_factory=list)
    total_weight: float = 0.0
    observations: int = 0

    def add(self, vector: Sequence[float], weight: float) -> None:
        if not vector or weight <= 0.0:
            return
        if not self.weighted_sum:
            self.weighted_sum = [0.0] * len(vector)
            self.weighted_sq_sum = [0.0] * len(vector)
        if len(vector) != len(self.weighted_sum):
            raise ValueError(
                f"inconsistent embedding dimension: expected {len(self.weighted_sum)}, got {len(vector)}"
            )
        for index, value in enumerate(vector):
            self.weighted_sum[index] += weight * value
            self.weighted_sq_sum[index] += weight * value * value
        self.total_weight += weight
        self.observations += 1

    def mean_std(self) -> List[float]:
        if self.total_weight <= 0.0:
            return []
        mean = [value / self.total_weight for value in self.weighted_sum]
        variance = [
            max(0.0, sq / self.total_weight - mu * mu)
            for sq, mu in zip(self.weighted_sq_sum, mean)
        ]
        return mean + [math.sqrt(value) for value in variance]

    def confidence(self, saturation: float = 4.0) -> float:
        if self.total_weight <= 0.0:
            return 0.0
        mean_weight = min(1.0, self.total_weight / max(1, self.observations))
        amount = 1.0 - math.exp(-self.observations / saturation)
        return max(0.0, min(1.0, mean_weight * amount))


@dataclass
class GeneProfile:
    values: List[float]
    confidence_sum: float = 0.0
    observations: int = 0

    def add(self, family: str, confidence: float) -> None:
        digest = hashlib.blake2b(family.encode("utf-8"), digest_size=16).digest()
        bucket = int.from_bytes(digest[:8], "little") % len(self.values)
        sign = 1.0 if digest[8] & 1 else -1.0
        self.values[bucket] += sign * confidence
        self.confidence_sum += confidence
        self.observations += 1

    def normalized(self) -> List[float]:
        norm = math.sqrt(sum(value * value for value in self.values))
        if norm <= 1e-12:
            return []
        return [value / norm for value in self.values]

    def confidence(self) -> float:
        if self.observations == 0:
            return 0.0
        quality = self.confidence_sum / self.observations
        amount = 1.0 - math.exp(-self.observations / 8.0)
        return max(0.0, min(1.0, quality * amount))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--dna-embeddings",
        type=Path,
        help="TSV: contig, embedding, optional window_id and confidence",
    )
    parser.add_argument(
        "--gene-hits",
        type=Path,
        help="TSV: contig, gene_id, family (or annotation), optional confidence",
    )
    parser.add_argument(
        "--protein-embeddings",
        type=Path,
        help="TSV: contig, protein_id, embedding, optional confidence; embedding is comma-separated",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        help="TSV: contig, taxonomy/lineage, optional confidence",
    )
    parser.add_argument("--gene-dim", type=int, default=256)
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
        if value is not None and value.strip():
            return value.strip()
    return default


def probability(raw: str, default: float = 1.0) -> float:
    if not raw or raw in {".", "NA", "na"}:
        return default
    value = float(raw)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"confidence must be within [0,1], got {raw!r}")
    return value


def vector(raw: str) -> List[float]:
    values = [float(value) for value in raw.split(",") if value.strip()]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("embedding must contain finite comma-separated floats")
    return values


def format_vector(values: Sequence[float]) -> str:
    return ",".join(f"{value:.7g}" for value in values)


def build(args: argparse.Namespace) -> Tuple[
    Dict[str, WeightedVectors],
    Dict[str, GeneProfile],
    Dict[str, WeightedVectors],
    Dict[str, Tuple[str, float]],
]:
    if args.gene_dim <= 0:
        raise ValueError("--gene-dim must be positive")

    dna: Dict[str, WeightedVectors] = defaultdict(WeightedVectors)
    genes: Dict[str, GeneProfile] = {}
    proteins: Dict[str, WeightedVectors] = defaultdict(WeightedVectors)
    taxonomy: Dict[str, Tuple[str, float]] = {}

    if args.dna_embeddings:
        for row in rows(args.dna_embeddings):
            contig = first(row, ("contig", "contig_id", "sequence"))
            embedding = first(row, ("embedding", "dna_embedding", "dna_lm_embedding", "dnabert_embedding"))
            if not contig or not embedding:
                continue
            confidence = probability(first(row, ("confidence", "score")), 1.0)
            dna[contig].add(vector(embedding), confidence)

    if args.gene_hits:
        for row in rows(args.gene_hits):
            contig = first(row, ("contig", "contig_id", "sequence"))
            family = first(row, ("family", "gene_family", "annotation", "function"))
            if not contig or not family or family in {".", "NA"}:
                continue
            confidence = probability(first(row, ("confidence", "score")), 1.0)
            profile = genes.setdefault(contig, GeneProfile([0.0] * args.gene_dim))
            profile.add(family, confidence)

    if args.protein_embeddings:
        for row in rows(args.protein_embeddings):
            contig = first(row, ("contig", "contig_id", "sequence"))
            embedding = first(row, ("embedding", "esm_embedding", "esmc_embedding"))
            if not contig or not embedding:
                continue
            confidence = probability(first(row, ("confidence", "score")), 1.0)
            proteins[contig].add(vector(embedding), confidence)

    if args.taxonomy:
        for row in rows(args.taxonomy):
            contig = first(row, ("contig", "contig_id", "sequence"))
            lineage = first(row, ("taxonomy", "lineage", "taxon"))
            if not contig or not lineage:
                continue
            confidence = probability(first(row, ("confidence", "taxonomy_confidence", "score")), 1.0)
            current = taxonomy.get(contig)
            if current is None or confidence > current[1]:
                taxonomy[contig] = (lineage, confidence)

    return dna, genes, proteins, taxonomy


def write_output(
    output: Path,
    dna: Dict[str, WeightedVectors],
    genes: Dict[str, GeneProfile],
    proteins: Dict[str, WeightedVectors],
    taxonomy: Dict[str, Tuple[str, float]],
) -> None:
    contigs = sorted(set(dna) | set(genes) | set(proteins) | set(taxonomy))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
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
        for contig in contigs:
            dna_feature = dna.get(contig)
            gene = genes.get(contig)
            protein = proteins.get(contig)
            lineage, tax_conf = taxonomy.get(contig, (".", 0.0))
            writer.writerow(
                [
                    contig,
                    lineage,
                    f"{tax_conf:.6f}",
                    format_vector(dna_feature.mean_std()) if dna_feature else ".",
                    f"{dna_feature.confidence(saturation=2.0):.6f}" if dna_feature else "0",
                    format_vector(gene.normalized()) if gene else ".",
                    f"{gene.confidence():.6f}" if gene else "0",
                    format_vector(protein.mean_std()) if protein else ".",
                    f"{protein.confidence():.6f}" if protein else "0",
                ]
            )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not any((args.dna_embeddings, args.gene_hits, args.protein_embeddings, args.taxonomy)):
        raise SystemExit(
            "provide at least one of --dna-embeddings, --gene-hits, --protein-embeddings, --taxonomy"
        )
    dna, genes, proteins, taxonomy = build(args)
    write_output(args.output, dna, genes, proteins, taxonomy)
    print(
        f"bridgebin-bio: contigs={len(set(dna) | set(genes) | set(proteins) | set(taxonomy))} "
        f"dna={len(dna)} gene={len(genes)} protein={len(proteins)} taxonomy={len(taxonomy)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
