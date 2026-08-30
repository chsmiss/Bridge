#!/usr/bin/env python3
"""Collapse exact duplicate protein sequences while preserving provenance.

PLASS may emit multiple FASTA records with exactly the same amino-acid sequence.
Those records are technical/support multiplicity, not independent homology alternatives.
If exact duplicates are counted independently, every informative amino-acid k-mer can be
mistaken for a highly repetitive k-mer and removed from the edge-evidence index.

This utility writes one representative record per exact sequence and an optional TSV that
maps representatives to all original record identifiers. It does not cluster non-identical
proteins and therefore does not erase genuine homolog ambiguity.
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
from typing import Iterator


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper().replace("*", "")
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"empty FASTA identifier in {path}")
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"sequence encountered before FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper().replace("*", "")


def write_record(handle, name: str, sequence: str, copies: int) -> None:
    handle.write(f">{name} exact_copies={copies}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mapping", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    sequence_to_ids: collections.OrderedDict[str, list[str]] = collections.OrderedDict()
    total = 0
    for protein_id, sequence in read_fasta(args.input):
        total += 1
        if not sequence:
            continue
        sequence_to_ids.setdefault(sequence, []).append(protein_id)

    if not sequence_to_ids:
        raise ValueError("protein FASTA contains no non-empty sequences")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    representatives: list[tuple[str, str, list[str]]] = []
    for index, (sequence, ids) in enumerate(sequence_to_ids.items(), start=1):
        representative = f"prot_cluster_{index:08d}"
        representatives.append((representative, sequence, ids))

    with args.output.open("w", encoding="utf-8") as handle:
        for representative, sequence, ids in representatives:
            write_record(handle, representative, sequence, len(ids))

    if args.mapping is not None:
        args.mapping.parent.mkdir(parents=True, exist_ok=True)
        with args.mapping.open("w", encoding="utf-8") as handle:
            handle.write("representative\texact_copies\toriginal_ids\n")
            for representative, _sequence, ids in representatives:
                handle.write(
                    f"{representative}\t{len(ids)}\t{','.join(ids)}\n"
                )

    copy_counts = [len(ids) for _representative, _sequence, ids in representatives]
    summary = {
        "input_records": total,
        "unique_exact_sequences": len(representatives),
        "collapsed_records": total - len(representatives),
        "largest_exact_copy_cluster": max(copy_counts),
        "mean_exact_copy_cluster": sum(copy_counts) / len(copy_counts),
    }
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
