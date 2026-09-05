#!/usr/bin/env python3
"""Encode prior contigs as exactly one synthetic physical fragment.

This is intentionally different from sliding virtual paired reads.  All accepted
prior contigs are concatenated into one R1 record, separated by N^target_k; R2
is a single N.  BridgeAsm deduplicates k-mers within one physical read pair
before incrementing fragment_count, so any target-k k-mer receives at most +1
synthetic fragment support globally.  N separators also break graph transitions
and read-thread paths between distinct prior contigs.

The intended use is multi-k carry-forward: the previous k supplies one unit of
structural prior, but a target-k k-mer must still occur in at least one physical
raw fragment to pass the production min-fragment-support=2 gate.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import low_abundance_rescue as lr


def build_prior_sequences(inputs: list[Path], min_length: int) -> tuple[list[str], dict[str, int]]:
    seen: set[str] = set()
    sequences: list[str] = []
    input_records = 0
    input_bases = 0
    short_records = 0
    duplicate_records = 0
    n_records = 0
    for path in inputs:
        if not path.exists() or path.stat().st_size == 0:
            continue
        for _name, seq0 in lr.fasta_records(path):
            input_records += 1
            seq = seq0.upper()
            input_bases += len(seq)
            if len(seq) < min_length:
                short_records += 1
                continue
            if "N" in seq:
                # A source-side N would make provenance of adjacent target-k
                # windows harder to reason about.  Keep this primitive strict.
                n_records += 1
                continue
            canonical = lr.canonical(seq)
            if canonical in seen:
                duplicate_records += 1
                continue
            seen.add(canonical)
            sequences.append(seq)
    return sequences, {
        "input_records": input_records,
        "input_bases": input_bases,
        "accepted_records": len(sequences),
        "accepted_bases": sum(map(len, sequences)),
        "short_records": short_records,
        "duplicate_records": duplicate_records,
        "n_records": n_records,
    }


def target_kmer_stats(sequences: list[str], k: int) -> dict[str, int]:
    observations = 0
    distinct: set[str] = set()
    repeated_within_prior = 0
    seen: set[str] = set()
    for seq in sequences:
        for mer in lr.kmers(seq, k):
            observations += 1
            if mer in seen:
                repeated_within_prior += 1
            else:
                seen.add(mer)
            distinct.add(mer)
    return {
        "target_kmer_observations": observations,
        "distinct_target_kmers": len(distinct),
        "repeated_target_kmer_observations": repeated_within_prior,
    }


def write_one_fragment(
    sequences: list[str],
    read1: Path,
    read2: Path,
    *,
    target_k: int,
    phred: int,
) -> dict[str, int]:
    if target_k <= 0:
        raise ValueError("target_k must be positive")
    if not 0 <= phred <= 40:
        raise ValueError("phred must be in 0..40")
    read1.parent.mkdir(parents=True, exist_ok=True)
    read2.parent.mkdir(parents=True, exist_ok=True)
    separator = "N" * target_k
    joined = separator.join(sequences)
    qchar = chr(33 + phred)
    with gzip.open(read1, "wt", compresslevel=3) as left, gzip.open(
        read2, "wt", compresslevel=3
    ) as right:
        if joined:
            left.write(f"@trusted_prior_k{target_k}/1\n{joined}\n+\n{qchar * len(joined)}\n")
            right.write(f"@trusted_prior_k{target_k}/2\nN\n+\n!\n")
    return {
        "synthetic_physical_fragments": int(bool(joined)),
        "synthetic_r1_bases": len(joined),
        "separator_bases": max(0, len(sequences) - 1) * target_k,
        "max_synthetic_fragment_support_per_kmer": int(bool(joined)),
        "phred": phred,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", type=Path)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--target-k", type=int, required=True)
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--phred", type=int, default=20)
    ap.add_argument("--stats-json", type=Path)
    args = ap.parse_args()

    sequences, source_stats = build_prior_sequences(args.inputs, args.min_length)
    stats = {
        "policy": "one_synthetic_physical_fragment_global",
        "target_k": args.target_k,
        "min_length": args.min_length,
        **source_stats,
        **target_kmer_stats(sequences, args.target_k),
        **write_one_fragment(
            sequences,
            args.read1,
            args.read2,
            target_k=args.target_k,
            phred=args.phred,
        ),
    }
    if args.stats_json:
        args.stats_json.parent.mkdir(parents=True, exist_ok=True)
        args.stats_json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
