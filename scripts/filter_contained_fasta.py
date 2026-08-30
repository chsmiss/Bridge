#!/usr/bin/env python3
"""Remove exact reverse-complement-aware contained FASTA records.

The filter is reference-free and conservative: a record is removed only after
its full sequence (or reverse complement) is found as an exact substring of an
already retained record of equal or greater length. A minimizer index is used
only to generate candidates; it can leave some redundant records behind, but
it cannot turn an approximate match into a removal.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Iterator

_BASE_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


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
    return sequence.translate(_COMPLEMENT)[::-1]


def canonical(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def canonical_kmer_values(sequence: str, k: int) -> list[int | None]:
    mask = (1 << (2 * k)) - 1
    top_shift = 2 * (k - 1)
    forward = 0
    reverse = 0
    valid = 0
    output: list[int | None] = []
    for base in sequence:
        bits = _BASE_BITS.get(base)
        if bits is None:
            forward = 0
            reverse = 0
            valid = 0
            output.append(None)
            continue
        forward = ((forward << 2) | bits) & mask
        reverse = (reverse >> 2) | ((3 - bits) << top_shift)
        valid += 1
        output.append(min(forward, reverse) if valid >= k else None)
    return output


def minimizers(sequence: str, k: int, window: int) -> list[int]:
    values = canonical_kmer_values(sequence, k)
    output: list[int] = []
    run: list[int] = []

    def emit_run(items: list[int]) -> None:
        if not items:
            return
        if len(items) <= window:
            value = min(items)
            if not output or output[-1] != value:
                output.append(value)
            return
        queue: deque[tuple[int, int]] = deque()
        for index, value in enumerate(items):
            while queue and queue[-1][1] >= value:
                queue.pop()
            queue.append((index, value))
            while queue[0][0] <= index - window:
                queue.popleft()
            if index >= window - 1:
                minimum = queue[0][1]
                if not output or output[-1] != minimum:
                    output.append(minimum)

    for value in values:
        if value is None:
            emit_run(run)
            run = []
        else:
            run.append(value)
    emit_run(run)
    return output


def write_fasta(path: Path, kept: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for index, (header, sequence) in enumerate(kept, start=1):
            handle.write(
                f">contained_filtered_{index:08d} len={len(sequence)} source={header}\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--min-length", type=int, default=0)
    parser.add_argument("--seed-k", type=int, default=21)
    parser.add_argument("--window", type=int, default=12)
    parser.add_argument("--candidate-minimizers", type=int, default=12)
    parser.add_argument("--removed-tsv", type=Path)
    parser.add_argument("--stats-json", type=Path)
    args = parser.parse_args()

    if args.seed_k <= 0 or args.window <= 0 or args.candidate_minimizers <= 0:
        raise ValueError("seed-k, window, and candidate-minimizers must be positive")

    unique: dict[str, str] = {}
    input_records = 0
    input_bases = 0
    for header, sequence in records(args.input):
        input_records += 1
        input_bases += len(sequence)
        if len(sequence) < args.min_length:
            continue
        if set(sequence) - set("ACGTN"):
            raise ValueError(f"unsupported FASTA character in {args.input}")
        unique.setdefault(canonical(sequence), header)

    ordered = sorted(
        ((header, sequence) for sequence, header in unique.items()),
        key=lambda record: (-len(record[1]), record[1]),
    )

    index: dict[int, list[int]] = defaultdict(list)
    kept: list[tuple[str, str]] = []
    removed: list[tuple[str, int, str, int]] = []

    for header, sequence in ordered:
        sequence_minimizers = list(
            dict.fromkeys(minimizers(sequence, args.seed_k, args.window))
        )
        ranked = sorted(
            (
                (len(index[value]), value)
                for value in sequence_minimizers
                if value in index
            ),
            key=lambda item: (item[0], item[1]),
        )
        candidates: list[int] = []
        seen_candidates: set[int] = set()
        for _frequency, value in ranked[: args.candidate_minimizers]:
            for candidate_id in index[value]:
                if candidate_id not in seen_candidates:
                    seen_candidates.add(candidate_id)
                    candidates.append(candidate_id)

        reverse: str | None = None
        container_id: int | None = None
        for candidate_id in candidates:
            target = kept[candidate_id][1]
            if sequence in target:
                container_id = candidate_id
                break
            if reverse is None:
                reverse = reverse_complement(sequence)
            if reverse in target:
                container_id = candidate_id
                break

        if container_id is not None:
            container_header, container_sequence = kept[container_id]
            removed.append(
                (header, len(sequence), container_header, len(container_sequence))
            )
            continue

        kept_id = len(kept)
        kept.append((header, sequence))
        for value in sequence_minimizers:
            index[value].append(kept_id)

    write_fasta(args.output, kept)

    if args.removed_tsv is not None:
        args.removed_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.removed_tsv.open("w") as handle:
            handle.write(
                "removed_header\tremoved_length\tcontainer_header\tcontainer_length\n"
            )
            for row in removed:
                handle.write("\t".join(map(str, row)) + "\n")

    stats = {
        "input_records": input_records,
        "input_bases": input_bases,
        "canonical_unique_records": len(unique),
        "kept_records": len(kept),
        "kept_bases": sum(len(sequence) for _header, sequence in kept),
        "contained_records_removed": len(removed),
        "seed_k": args.seed_k,
        "window": args.window,
        "candidate_minimizers": args.candidate_minimizers,
    }
    if args.stats_json is not None:
        args.stats_json.parent.mkdir(parents=True, exist_ok=True)
        args.stats_json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
