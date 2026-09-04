#!/usr/bin/env python3
"""Extend an immutable seed assembly with conservative exact end overlaps.

The seed contigs are never replaced, split, or deleted.  A higher-k candidate is
allowed to extend one seed end only when the exact suffix-prefix overlap is the
unique best choice in both directions.  Candidates that reciprocally match more
than one physical seed end are rejected: those are potential bridges/repeat
resolutions and require stronger physical-fragment evidence than exact sequence
identity alone.

This is intentionally conservative.  It is designed for the Stage24 -> k31
carry-forward path where preserving rescued k21 sequence and avoiding cross-k
rearrangements matter more than maximizing raw N50.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterator

TRANS = bytes.maketrans(b"ACGTN", b"TGCAN")


def records(path: Path) -> Iterator[tuple[str, bytes]]:
    header: str | None = None
    chunks: list[bytes] = []
    with path.open("rb") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if header is not None:
                    yield header, b"".join(chunks).upper()
                header = line[1:].decode("utf-8", "replace")
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"sequence before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        yield header, b"".join(chunks).upper()


def rc(sequence: bytes) -> bytes:
    return sequence.translate(TRANS)[::-1]


@dataclass(frozen=True)
class ExtensionMatch:
    seed_state: int
    candidate_state: int
    overlap: int
    extension: int


def top_unique(values: dict[int, int], margin: int) -> tuple[int, int] | None:
    if not values:
        return None
    ranked = sorted(values.items(), key=lambda item: (-item[1], item[0]))
    if len(ranked) > 1:
        best = ranked[0][1]
        second = ranked[1][1]
        if best == second or best - second < margin:
            return None
    return ranked[0]


def n50(sequences: list[bytes]) -> int:
    lengths = sorted((len(sequence) for sequence in sequences), reverse=True)
    if not lengths:
        return 0
    threshold = (sum(lengths) + 1) // 2
    total = 0
    for length in lengths:
        total += length
        if total >= threshold:
            return length
    return 0


def find_reciprocal_extensions(
    seeds: list[bytes],
    candidates: list[bytes],
    *,
    min_overlap: int,
    overlap_margin: int,
    seed_length: int,
    min_extension: int,
    max_seed_occurrences: int,
) -> tuple[list[ExtensionMatch], int]:
    oriented_seeds: list[bytes] = []
    oriented_candidates: list[bytes] = []
    for sequence in seeds:
        oriented_seeds.extend((sequence, rc(sequence)))
    for sequence in candidates:
        oriented_candidates.extend((sequence, rc(sequence)))

    prefix: dict[bytes, list[int]] = defaultdict(list)
    for state, sequence in enumerate(oriented_candidates):
        if len(sequence) < seed_length:
            continue
        key = sequence[:seed_length]
        if set(key) <= set(b"ACGT"):
            prefix[key].append(state)
    prefix = {
        key: states
        for key, states in prefix.items()
        if len(states) <= max_seed_occurrences
    }

    outgoing: list[dict[int, int]] = [dict() for _ in oriented_seeds]
    incoming: list[dict[int, int]] = [dict() for _ in oriented_candidates]
    checks = 0

    for seed_state, seed in enumerate(oriented_seeds):
        max_start = len(seed) - min_overlap
        if max_start < 0:
            continue
        for start in range(max_start + 1):
            key = seed[start : start + seed_length]
            states = prefix.get(key)
            if not states:
                continue
            overlap = len(seed) - start
            suffix = seed[start:]
            for candidate_state in states:
                candidate = oriented_candidates[candidate_state]
                extension = len(candidate) - overlap
                if extension < min_extension:
                    continue
                checks += 1
                if not candidate.startswith(suffix):
                    continue
                old = outgoing[seed_state].get(candidate_state, 0)
                if overlap > old:
                    outgoing[seed_state][candidate_state] = overlap
                old = incoming[candidate_state].get(seed_state, 0)
                if overlap > old:
                    incoming[candidate_state][seed_state] = overlap

    best_out = [top_unique(values, overlap_margin) for values in outgoing]
    best_in = [top_unique(values, overlap_margin) for values in incoming]
    reciprocal: list[ExtensionMatch] = []
    for seed_state, item in enumerate(best_out):
        if item is None:
            continue
        candidate_state, overlap = item
        if best_in[candidate_state] != (seed_state, overlap):
            continue
        extension = len(oriented_candidates[candidate_state]) - overlap
        reciprocal.append(
            ExtensionMatch(seed_state, candidate_state, overlap, extension)
        )
    return reciprocal, checks


def seed_locked_extensions(
    seed_records: list[tuple[str, bytes]],
    candidate_records: list[tuple[str, bytes]],
    *,
    min_overlap: int,
    overlap_margin: int,
    seed_length: int,
    min_extension: int,
    max_seed_occurrences: int,
) -> tuple[list[tuple[str, bytes]], dict[str, int | float]]:
    seeds = [sequence for _header, sequence in seed_records]
    candidates = [sequence for _header, sequence in candidate_records]
    reciprocal, checks = find_reciprocal_extensions(
        seeds,
        candidates,
        min_overlap=min_overlap,
        overlap_margin=overlap_margin,
        seed_length=seed_length,
        min_extension=min_extension,
        max_seed_occurrences=max_seed_occurrences,
    )

    # A physical candidate with reciprocal matches to multiple seed ends is a
    # bridge/repeat-resolution hypothesis.  Exact overlap alone is insufficient
    # evidence for that operation, so reject all of its matches here.
    physical_candidate_matches = Counter(
        match.candidate_state // 2 for match in reciprocal
    )
    accepted = [
        match
        for match in reciprocal
        if physical_candidate_matches[match.candidate_state // 2] == 1
    ]

    left: dict[int, bytes] = {}
    right: dict[int, bytes] = {}
    for match in accepted:
        seed_id = match.seed_state // 2
        candidate_id = match.candidate_state // 2
        candidate = candidates[candidate_id]
        oriented_candidate = (
            candidate if match.candidate_state % 2 == 0 else rc(candidate)
        )
        extra = oriented_candidate[match.overlap :]
        if match.seed_state % 2 == 0:
            right[seed_id] = extra
        else:
            left[seed_id] = rc(extra)

    output: list[tuple[str, bytes]] = []
    for seed_id, (header, sequence) in enumerate(seed_records):
        left_extra = left.get(seed_id, b"")
        right_extra = right.get(seed_id, b"")
        merged = left_extra + sequence + right_extra
        output.append(
            (
                f"seed_locked_{seed_id + 1:08d} len={len(merged)} "
                f"left_extension={len(left_extra)} right_extension={len(right_extra)} "
                f"source={header}",
                merged,
            )
        )

    output_sequences = [sequence for _header, sequence in output]
    seed_bp = sum(len(sequence) for sequence in seeds)
    output_bp = sum(len(sequence) for sequence in output_sequences)
    stats: dict[str, int | float] = {
        "seed_records": len(seeds),
        "candidate_records": len(candidates),
        "candidate_checks": checks,
        "reciprocal_matches": len(reciprocal),
        "ambiguous_bridge_candidates_rejected": sum(
            count > 1 for count in physical_candidate_matches.values()
        ),
        "accepted_extensions": len(accepted),
        "left_extensions": len(left),
        "right_extensions": len(right),
        "seed_bp": seed_bp,
        "output_bp": output_bp,
        "added_bp": output_bp - seed_bp,
        "seed_n50": n50(seeds),
        "output_n50": n50(output_sequences),
    }
    return output, stats


def write_fasta(path: Path, output: list[tuple[str, bytes]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        for header, sequence in output:
            handle.write(b">" + header.encode() + b"\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + b"\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("seed", type=Path)
    parser.add_argument("candidates", type=Path)
    parser.add_argument("--min-overlap", type=int, default=200)
    parser.add_argument("--overlap-margin", type=int, default=30)
    parser.add_argument("--seed-length", type=int, default=31)
    parser.add_argument("--min-extension", type=int, default=20)
    parser.add_argument("--max-seed-occurrences", type=int, default=64)
    parser.add_argument("--stats-json", type=Path)
    args = parser.parse_args()

    if args.seed_length > args.min_overlap:
        raise SystemExit("--seed-length cannot exceed --min-overlap")
    if args.min_overlap < 1 or args.min_extension < 1:
        raise SystemExit("overlap and extension thresholds must be positive")

    seed_records = list(records(args.seed))
    candidate_records = list(records(args.candidates))
    if not seed_records:
        raise SystemExit("seed FASTA is empty")
    if not candidate_records:
        raise SystemExit("candidate FASTA is empty")

    output, stats = seed_locked_extensions(
        seed_records,
        candidate_records,
        min_overlap=args.min_overlap,
        overlap_margin=args.overlap_margin,
        seed_length=args.seed_length,
        min_extension=args.min_extension,
        max_seed_occurrences=args.max_seed_occurrences,
    )
    write_fasta(args.output, output)
    if args.stats_json:
        args.stats_json.parent.mkdir(parents=True, exist_ok=True)
        args.stats_json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
