#!/usr/bin/env python3
"""Emit simple canonical k-mer frequency embeddings for Biological Brain benchmarks.

This is intentionally boring. It provides a fair sequence-composition baseline for
asking whether a DNA foundation model adds genome-identity information beyond what the
existing binner can obtain cheaply from short k-mers.

Output columns are ``contig``, ``kmer_embedding``, ``kmer_confidence``, ``k``.
Reverse complements share one canonical bucket; the vector is L2 normalized.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("-k", type=int, default=5)
    return parser.parse_args(argv)


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"{path}: sequence before FASTA header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def encode(kmer: str) -> Optional[int]:
    value = 0
    for base in kmer:
        digit = BASE.get(base)
        if digit is None:
            return None
        value = (value << 2) | digit
    return value


def reverse_complement_code(kmer: str) -> Optional[int]:
    value = 0
    for base in reversed(kmer):
        digit = BASE.get(base)
        if digit is None:
            return None
        value = (value << 2) | (3 - digit)
    return value


def vector(sequence: str, k: int) -> Tuple[List[float], int]:
    size = 4**k
    counts = [0.0] * size
    observed = 0
    for start in range(0, max(0, len(sequence) - k + 1)):
        kmer = sequence[start : start + k]
        forward = encode(kmer)
        reverse = reverse_complement_code(kmer)
        if forward is None or reverse is None:
            continue
        counts[min(forward, reverse)] += 1.0
        observed += 1
    norm = math.sqrt(sum(value * value for value in counts))
    if norm > 0.0:
        counts = [value / norm for value in counts]
    return counts, observed


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not 2 <= args.k <= 7:
        raise SystemExit("-k must be between 2 and 7 for this diagnostic")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "kmer_embedding", "kmer_confidence", "k"])
        for contig, sequence in read_fasta(args.contigs):
            embedding, observed = vector(sequence, args.k)
            confidence = 1.0 - math.exp(-observed / max(1.0, 20.0 * (4**args.k)))
            writer.writerow(
                [
                    contig,
                    ",".join(f"{value:.7g}" for value in embedding),
                    f"{confidence:.6f}",
                    args.k,
                ]
            )
            written += 1
    print(f"bridgebin-kmer-baseline: contigs={written} k={args.k} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
