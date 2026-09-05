#!/usr/bin/env python3
"""Truth evaluator for synthetic major/minor strain assemblies.

Primary genome fractions are measured from the primary contig FASTA only.
Minor-allele recall is measured across primary contigs plus optional variant or
haplotig FASTAs, so retaining a true alternate allele is not confused with
forcing it into the primary backbone.
"""
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


def canonical(sequence: str) -> str:
    reverse = revcomp(sequence)
    return min(sequence, reverse)


def covered(reference: str, contigs: list[str], seed: int = 21) -> list[bool]:
    """Conservative exact-match breadth proxy for deterministic simulations."""
    mask = [False] * len(reference)
    seed_index: dict[str, list[int]] = {}
    if len(reference) < seed:
        return mask
    for position in range(len(reference) - seed + 1):
        kmer = reference[position : position + seed]
        seed_index.setdefault(kmer, []).append(position)

    for contig in contigs:
        for orientation in (contig, revcomp(contig)):
            if len(orientation) < seed:
                continue
            checked: set[int] = set()
            for offset in range(len(orientation) - seed + 1):
                kmer = orientation[offset : offset + seed]
                for reference_position in seed_index.get(kmer, []):
                    start = reference_position - offset
                    if start in checked or start < 0:
                        continue
                    checked.add(start)
                    end = min(len(reference), start + len(orientation))
                    observed = orientation[: end - start]
                    if reference[start:end] == observed:
                        mask[start:end] = [True] * (end - start)
    return mask


def output_kmers(sequences: list[str], k: int) -> set[str]:
    kmers: set[str] = set()
    for sequence in sequences:
        if len(sequence) < k:
            continue
        for position in range(len(sequence) - k + 1):
            window = sequence[position : position + k]
            if "N" not in window:
                kmers.add(canonical(window))
    return kmers


def allele_window(sequence: str, position: int, k: int) -> str | None:
    if len(sequence) < k or position < 0 or position >= len(sequence):
        return None
    left = max(0, position - k // 2)
    left = min(left, len(sequence) - k)
    window = sequence[left : left + k]
    return window if len(window) == k else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--major", type=Path, required=True)
    parser.add_argument("--minor", type=Path, required=True)
    parser.add_argument(
        "--contigs",
        type=Path,
        required=True,
        help="Primary contig FASTA used for genome-fraction estimates",
    )
    parser.add_argument(
        "--additional",
        type=Path,
        action="append",
        default=[],
        help="Additional variant/haplotig FASTA used only for allele recall",
    )
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--allele-k", type=int, default=31)
    args = parser.parse_args()

    major = read_fasta(args.major)[0]
    minor = read_fasta(args.minor)[0]
    primary = read_fasta(args.contigs)
    additional = [sequence for path in args.additional for sequence in read_fasta(path)]
    all_sequences = primary + additional

    major_mask = covered(major, primary)
    minor_mask = covered(minor, primary)
    variant_positions = [
        int(line)
        for line in args.variants.read_text().splitlines()[1:]
        if line.strip()
    ]

    allele_k = min(args.allele_k, len(minor))
    if allele_k % 2 == 0:
        allele_k -= 1
    observed_kmers = output_kmers(all_sequences, allele_k)
    informative = 0
    recalled = 0
    for position in variant_positions:
        minor_window = allele_window(minor, position, allele_k)
        major_window = allele_window(major, position, allele_k)
        if minor_window is None or major_window is None or minor_window == major_window:
            continue
        informative += 1
        if canonical(minor_window) in observed_kmers:
            recalled += 1

    print(f"primary_major_fraction\t{sum(major_mask) / len(major):.6f}")
    print(f"primary_minor_fraction\t{sum(minor_mask) / len(minor):.6f}")
    print(f"minor_variant_recall\t{recalled / max(1, informative):.6f}")
    print(f"minor_variants_informative\t{informative}")
    print(f"minor_variants_recalled\t{recalled}")
    print(f"primary_contigs\t{len(primary)}")
    print(f"primary_assembled_bp\t{sum(map(len, primary))}")
    print(f"additional_sequences\t{len(additional)}")
    print(f"all_output_bp\t{sum(map(len, all_sequences))}")


if __name__ == "__main__":
    main()
