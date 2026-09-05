#!/usr/bin/env python3
"""Deterministic paired-end metagenome simulator for BridgeAsm smoke benchmarks."""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

DNA = "ACGT"


def revcomp(sequence: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


def mutate(sequence: str, rate: float, rng: random.Random) -> tuple[str, list[int]]:
    output = list(sequence)
    positions: list[int] = []
    for index, base in enumerate(output):
        if rng.random() < rate:
            output[index] = rng.choice([candidate for candidate in DNA if candidate != base])
            positions.append(index)
    return "".join(output), positions


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return path.open("w")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome-length", type=int, default=50_000)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--insert-mean", type=int, default=350)
    parser.add_argument("--major-depth", type=float, default=30.0)
    parser.add_argument("--minor-depth", type=float, default=5.0)
    parser.add_argument("--strain-divergence", type=float, default=0.005)
    parser.add_argument("--error-rate", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=43)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    major = "".join(rng.choice(DNA) for _ in range(args.genome_length))
    minor, variants = mutate(major, args.strain_divergence, rng)
    (args.output / "major.fasta").write_text(f">major\n{major}\n")
    (args.output / "minor.fasta").write_text(f">minor\n{minor}\n")
    (args.output / "minor_variants.tsv").write_text(
        "position\n" + "".join(f"{position}\n" for position in variants)
    )

    r1_path = args.output / "reads_R1.fastq.gz"
    r2_path = args.output / "reads_R2.fastq.gz"
    with open_text(r1_path) as r1, open_text(r2_path) as r2:
        read_id = 0
        for label, genome, depth in [
            ("major", major, args.major_depth),
            ("minor", minor, args.minor_depth),
        ]:
            pairs = int(depth * len(genome) / (2 * args.read_length))
            for _ in range(pairs):
                insert = max(2 * args.read_length, int(rng.gauss(args.insert_mean, 25)))
                if insert >= len(genome):
                    continue
                start = rng.randrange(0, len(genome) - insert)
                left = list(genome[start : start + args.read_length])
                right_start = start + insert - args.read_length
                right = list(revcomp(genome[right_start : right_start + args.read_length]))
                for read in (left, right):
                    for index, base in enumerate(read):
                        if rng.random() < args.error_rate:
                            read[index] = rng.choice(
                                [candidate for candidate in DNA if candidate != base]
                            )
                quality = "I" * args.read_length
                r1.write(f"@{label}_{read_id}/1\n{''.join(left)}\n+\n{quality}\n")
                r2.write(f"@{label}_{read_id}/2\n{''.join(right)}\n+\n{quality}\n")
                read_id += 1


if __name__ == "__main__":
    main()
