#!/usr/bin/env python3
"""Collapse technical duplicate/contained PLASS protein records with provenance.

PLASS can emit many records that are exact copies or exact substrings of a longer assembled
protein from the same read set. Counting those records as independent alternatives inflates
k-mer occurrence and homology ambiguity even though they encode the same collinear path.

This utility always collapses exact duplicate sequences. With --collapse-contained it also
collapses a sequence only when the *entire amino-acid sequence* is an exact substring of an
already retained longer sequence. Non-identical, non-contained homologs remain separate, so
true paralog/repeat ambiguity is preserved.
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
from typing import Iterator


CONTAINMENT_SEED = 12


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


def write_record(handle, name: str, sequence: str, copies: int, members: int) -> None:
    handle.write(f">{name} exact_copies={copies} member_records={members}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def collapse_contained(
    exact_records: list[tuple[str, list[str]]],
) -> tuple[list[tuple[str, list[str]]], int]:
    """Merge exact full-sequence substrings into longer retained proteins.

    Records are processed longest first. A 12-aa seed index limits substring checks to
    plausible containers rather than comparing every pair. Sequences shorter than the seed
    are retained because collapsing very short peptides would be biologically unsafe.
    """

    ordered = sorted(exact_records, key=lambda item: (-len(item[0]), item[0]))
    retained: list[list[object]] = []  # [sequence, original_ids]
    seed_to_retained: dict[str, list[int]] = collections.defaultdict(list)
    collapsed = 0

    for sequence, original_ids in ordered:
        container_index: int | None = None
        if len(sequence) >= CONTAINMENT_SEED:
            seed = sequence[:CONTAINMENT_SEED]
            for index in seed_to_retained.get(seed, []):
                container = retained[index][0]
                assert isinstance(container, str)
                if sequence in container:
                    container_index = index
                    break

        if container_index is not None:
            members = retained[container_index][1]
            assert isinstance(members, list)
            members.extend(original_ids)
            collapsed += 1
            continue

        retained_index = len(retained)
        retained.append([sequence, list(original_ids)])
        if len(sequence) >= CONTAINMENT_SEED:
            seen_seeds = {
                sequence[start : start + CONTAINMENT_SEED]
                for start in range(len(sequence) - CONTAINMENT_SEED + 1)
            }
            for seed in seen_seeds:
                seed_to_retained[seed].append(retained_index)

    return [(str(sequence), list(ids)) for sequence, ids in retained], collapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mapping", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument(
        "--collapse-contained",
        action="store_true",
        help="collapse exact full-sequence substrings into retained longer proteins",
    )
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

    exact_records = list(sequence_to_ids.items())
    unique_exact_sequences = len(exact_records)
    contained_sequences_collapsed = 0
    if args.collapse_contained:
        exact_records, contained_sequences_collapsed = collapse_contained(exact_records)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    representatives: list[tuple[str, str, list[str]]] = []
    for index, (sequence, ids) in enumerate(exact_records, start=1):
        representative = f"prot_cluster_{index:08d}"
        representatives.append((representative, sequence, ids))

    with args.output.open("w", encoding="utf-8") as handle:
        for representative, sequence, ids in representatives:
            write_record(handle, representative, sequence, len(ids), len(ids))

    if args.mapping is not None:
        args.mapping.parent.mkdir(parents=True, exist_ok=True)
        with args.mapping.open("w", encoding="utf-8") as handle:
            handle.write("representative\tmember_records\toriginal_ids\n")
            for representative, _sequence, ids in representatives:
                handle.write(f"{representative}\t{len(ids)}\t{','.join(ids)}\n")

    member_counts = [len(ids) for _representative, _sequence, ids in representatives]
    summary = {
        "input_records": total,
        "unique_exact_sequences": unique_exact_sequences,
        "exact_duplicate_records_collapsed": total - unique_exact_sequences,
        "contained_sequences_collapsed": contained_sequences_collapsed,
        "retained_sequences": len(representatives),
        "largest_member_cluster": max(member_counts),
        "mean_member_cluster": sum(member_counts) / len(member_counts),
        "containment_enabled": args.collapse_contained,
        "containment_seed_aa": CONTAINMENT_SEED,
    }
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
