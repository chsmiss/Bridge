#!/usr/bin/env python3
"""Merge FASTA files with exact reverse-complement-aware deduplication.

This utility deliberately performs only exact sequence deduplication. It does
not use a reference, align overlapping contigs, or create new joins. It is used
to measure how much complementary sequence independent k values recover.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator


def records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"sequence before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def canonical(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--min-length", type=int, default=0)
    args = parser.parse_args()

    sources: dict[str, set[str]] = {}
    for path in args.inputs:
        for _header, sequence in records(path):
            if len(sequence) < args.min_length:
                continue
            if set(sequence) - set("ACGTN"):
                raise ValueError(f"unsupported FASTA character in {path}")
            sources.setdefault(canonical(sequence), set()).add(path.name)

    ordered = sorted(sources, key=lambda sequence: (-len(sequence), sequence))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for index, sequence in enumerate(ordered, start=1):
            source_names = ",".join(sorted(sources[sequence]))
            handle.write(
                f">merged_{index:08d} len={len(sequence)} sources={source_names}\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    print(
        f"merged {len(args.inputs)} FASTA files into {len(ordered)} exact canonical records"
    )


if __name__ == "__main__":
    main()
