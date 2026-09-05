#!/usr/bin/env python3
"""Convert FASTA contigs into a link-free GFA backbone.

Each FASTA record becomes one GFA segment. Existing contigs are therefore never
split by the protein-guided path-cover stage; only new, explicitly supported
joins can be added between them.
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

COVERAGE_RE = re.compile(
    r"(?:^|\s)(?:cov|coverage|depth)=([0-9]+(?:\.[0-9]+)?)",
    re.IGNORECASE,
)
IUPAC_DNA = frozenset("ACGTNRYKMSWBDHV")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def fasta_records(path: Path) -> Iterator[tuple[str, str, str]]:
    name: str | None = None
    header = ""
    sequence: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(sequence).upper()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty FASTA header at line {line_number}")
                name = header.split()[0]
                sequence = []
            else:
                if name is None:
                    raise ValueError(
                        f"sequence before first FASTA header at line {line_number}"
                    )
                sequence.append(line)
    if name is not None:
        yield name, header, "".join(sequence).upper()


def coverage_from_header(header: str, default: float) -> float:
    match = COVERAGE_RE.search(header)
    if match is None:
        return default
    value = float(match.group(1))
    if value <= 0.0:
        return default
    return value


def validate_sequence(name: str, sequence: str) -> None:
    if not sequence:
        raise ValueError(f"FASTA record {name!r} has an empty sequence")
    invalid = sorted(set(sequence) - IUPAC_DNA)
    if invalid:
        preview = "".join(invalid[:12])
        raise ValueError(f"FASTA record {name!r} contains invalid DNA symbols: {preview}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path)
    parser.add_argument("gfa", type=Path)
    parser.add_argument("--default-coverage", type=float, default=1.0)
    parser.add_argument("--min-length", type=int, default=1)
    args = parser.parse_args()

    if args.default_coverage <= 0.0:
        parser.error("--default-coverage must be positive")
    if args.min_length < 1:
        parser.error("--min-length must be at least one")

    args.gfa.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    written = 0
    total_bp = 0
    skipped = 0
    with args.gfa.open("w") as writer:
        writer.write("H\tVN:Z:1.0\n")
        for name, header, sequence in fasta_records(args.fasta):
            if name in seen:
                raise ValueError(f"duplicate FASTA identifier: {name}")
            seen.add(name)
            validate_sequence(name, sequence)
            if len(sequence) < args.min_length:
                skipped += 1
                continue
            coverage = coverage_from_header(header, args.default_coverage)
            writer.write(
                f"S\t{name}\t{sequence}\tLN:i:{len(sequence)}\tKC:f:{coverage:.6f}\n"
            )
            written += 1
            total_bp += len(sequence)

    if written == 0:
        raise ValueError("no FASTA records passed --min-length")
    print(f"segments\t{written}")
    print(f"total_bp\t{total_bp}")
    print(f"skipped\t{skipped}")
    print(f"output\t{args.gfa}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
