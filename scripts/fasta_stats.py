#!/usr/bin/env python3
"""Report deterministic FASTA length statistics as TSV and JSON."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def lengths(path: Path) -> list[int]:
    output: list[int] = []
    current = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    output.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        output.append(current)
    return output


def n50(values: list[int]) -> int:
    if not values:
        return 0
    total = sum(values)
    cumulative = 0
    for value in sorted(values, reverse=True):
        cumulative += value
        if cumulative * 2 >= total:
            return value
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    values = lengths(args.fasta)
    stats = {
        "file": str(args.fasta),
        "sequences": len(values),
        "total_bp": sum(values),
        "n50": n50(values),
        "largest": max(values, default=0),
        "ge_500": sum(value >= 500 for value in values),
        "ge_1000": sum(value >= 1000 for value in values),
        "ge_5000": sum(value >= 5000 for value in values),
    }
    for key, value in stats.items():
        print(f"{key}\t{value}")
    if args.json:
        args.json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
