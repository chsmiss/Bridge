#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


def fasta_lengths(path: Path) -> List[int]:
    lengths: List[int] = []
    current = 0
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
    if current:
        lengths.append(current)
    return lengths


def gfa_lengths(path: Path) -> List[int]:
    lengths: List[int] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == "S" and fields[2] != "*":
                lengths.append(len(fields[2]))
    return lengths


def n_stat(lengths: Iterable[int], fraction: float) -> int:
    ordered = sorted((length for length in lengths if length > 0), reverse=True)
    total = sum(ordered)
    threshold = total * fraction
    cumulative = 0
    for length in ordered:
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def summarize(lengths: List[int], minimum: int) -> Dict[str, int]:
    kept = [length for length in lengths if length >= minimum]
    return {
        "sequences": len(kept),
        "total_bp": sum(kept),
        "largest": max(kept, default=0),
        "n50": n_stat(kept, 0.50),
        "n90": n_stat(kept, 0.90),
        "minimum": minimum,
    }


def evidence_summary(path: Optional[Path]) -> Dict[str, int]:
    if path is None:
        return {}
    counts: Dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            label = row.get("breakpoint_class", "unknown") or "unknown"
            counts[label] = counts.get(label, 0) + 1
    return counts


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--format", choices=("fasta", "gfa"), required=True)
    parser.add_argument("--minimum", type=int, default=200)
    parser.add_argument("--evidence", type=Path)
    parser.add_argument("--label", default="assembly")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    lengths = fasta_lengths(args.input) if args.format == "fasta" else gfa_lengths(args.input)
    result = {
        "label": args.label,
        **summarize(lengths, args.minimum),
        "evidence_classes": evidence_summary(args.evidence),
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
