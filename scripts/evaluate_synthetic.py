#!/usr/bin/env python3
"""Small truth evaluator: reference breadth and minor-allele presence."""
from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> list[str]:
    records: list[str] = []
    current: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current:
                records.append("".join(current))
                current = []
        else:
            current.append(line.strip().upper())
    if current:
        records.append("".join(current))
    return records


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def covered(reference: str, contigs: list[str], seed: int = 21) -> list[bool]:
    mask = [False] * len(reference)
    index: dict[str, list[int]] = {}
    for position in range(0, len(reference) - seed + 1):
        kmer = reference[position : position + seed]
        index.setdefault(kmer, []).append(position)
        index.setdefault(revcomp(kmer), []).append(position)
    for contig in contigs:
        for orientation in (contig, revcomp(contig)):
            for offset in range(0, max(0, len(orientation) - seed + 1)):
                kmer = orientation[offset : offset + seed]
                for position in index.get(kmer, []):
                    start = max(0, position - offset)
                    end = min(len(reference), start + len(orientation))
                    if reference[start:end] == orientation[: end - start]:
                        for base in range(start, end):
                            mask[base] = True
    return mask


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--major", type=Path, required=True)
    parser.add_argument("--minor", type=Path, required=True)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--variants", type=Path, required=True)
    args = parser.parse_args()

    major = read_fasta(args.major)[0]
    minor = read_fasta(args.minor)[0]
    contigs = read_fasta(args.contigs)
    major_mask = covered(major, contigs)
    minor_mask = covered(minor, contigs)
    variants = [int(line) for line in args.variants.read_text().splitlines()[1:] if line]
    minor_recalled = sum(1 for position in variants if minor_mask[position])
    print(f"major_fraction\t{sum(major_mask) / len(major):.6f}")
    print(f"minor_fraction\t{sum(minor_mask) / len(minor):.6f}")
    print(f"minor_variant_recall\t{minor_recalled / max(1, len(variants)):.6f}")
    print(f"contigs\t{len(contigs)}")
    print(f"assembled_bp\t{sum(map(len, contigs))}")


if __name__ == "__main__":
    main()
