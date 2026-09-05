#!/usr/bin/env python3
"""Segment-level replacement merge for graph backbone plus recovery contigs.

The new graph backbone is authoritative where it already represents recovery
sequence. Recovery contigs that are mostly represented are *not* discarded as
whole records: only their novel intervals are retained, with a small represented
flank so exact graph/sequence overlap can still be recognized downstream.

This keeps the duplication benefit of backbone replacement without deleting a
short unique extension, bubble allele, or low-coverage island merely because it
sits on an otherwise redundant recovery contig.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

BASE = {65: 0, 67: 1, 71: 2, 84: 3, 97: 0, 99: 1, 103: 2, 116: 3}


def fasta(path: Path):
    name = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def canonical_keys_with_pos(seq: str, k: int):
    mask = (1 << (2 * k)) - 1
    forward = reverse = valid = 0
    shift = 2 * (k - 1)
    for i, byte in enumerate(seq.encode("ascii", "ignore")):
        value = BASE.get(byte)
        if value is None:
            forward = reverse = valid = 0
            continue
        forward = ((forward << 2) | value) & mask
        reverse = (reverse >> 2) | ((3 - value) << shift)
        valid += 1
        if valid >= k:
            yield i - k + 1, min(forward, reverse)


def canonical_keys(seq: str, k: int):
    for _, key in canonical_keys_with_pos(seq, k):
        yield key


def merge_intervals(intervals: list[tuple[int, int]], gap: int = 0) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        old_start, old_end = merged[-1]
        if start <= old_end + gap:
            merged[-1] = (old_start, max(old_end, end))
        else:
            merged.append((start, end))
    return merged


def novel_intervals(
    seq: str,
    backbone_keys: set[int],
    k: int,
    min_novel_kmers: int,
    join_gap_kmers: int,
    flank: int,
) -> tuple[list[tuple[int, int]], int, int]:
    positions = list(canonical_keys_with_pos(seq, k))
    if not positions:
        return [], 0, 0
    novel_starts = [pos for pos, key in positions if key not in backbone_keys]
    represented = len(positions) - len(novel_starts)
    if not novel_starts:
        return [], len(positions), represented

    runs: list[tuple[int, int, int]] = []
    run_start = run_prev = novel_starts[0]
    count = 1
    for pos in novel_starts[1:]:
        if pos <= run_prev + join_gap_kmers + 1:
            run_prev = pos
            count += 1
            continue
        runs.append((run_start, run_prev + k, count))
        run_start = run_prev = pos
        count = 1
    runs.append((run_start, run_prev + k, count))

    intervals = []
    for start, end, count in runs:
        if count < min_novel_kmers:
            continue
        intervals.append((max(0, start - flank), min(len(seq), end + flank)))
    return merge_intervals(intervals, gap=2 * flank), len(positions), represented


def write_record(handle, name: str, seq: str) -> None:
    handle.write(f">{name} len={len(seq)}\n")
    for start in range(0, len(seq), 80):
        handle.write(seq[start : start + 80] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", type=Path, required=True)
    ap.add_argument("--recovery", type=Path, required=True)
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--stats-json", type=Path)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("--segment-fraction", type=float, default=0.85,
                    help="segment a recovery contig when this fraction of its k-mers is represented")
    ap.add_argument("--min-informative-kmers", type=int, default=20)
    ap.add_argument("--min-novel-kmers", type=int, default=4)
    ap.add_argument("--join-gap-kmers", type=int, default=2)
    ap.add_argument("--flank", type=int, default=30)
    ap.add_argument("--min-segment-length", type=int, default=80)
    args = ap.parse_args()

    backbone = list(fasta(args.backbone))
    recovery = list(fasta(args.recovery))
    backbone_keys: set[int] = set()
    for _, seq in backbone:
        backbone_keys.update(canonical_keys(seq, args.k))

    rows: list[tuple[object, ...]] = []
    kept_whole: list[tuple[str, str]] = []
    kept_segments: list[tuple[str, str]] = []
    fully_redundant_records = 0
    segmented_records = 0
    segmented_source_bases = 0
    novel_segment_bases = 0

    for name, seq in recovery:
        positions = list(canonical_keys_with_pos(seq, args.k))
        informative = len(positions)
        represented = sum(key in backbone_keys for _, key in positions)
        fraction = represented / max(1, informative)
        if informative < args.min_informative_kmers or fraction < args.segment_fraction:
            kept_whole.append((name, seq))
            rows.append((name, len(seq), informative, represented, fraction, "whole", 0, len(seq)))
            continue

        intervals, _, _ = novel_intervals(
            seq,
            backbone_keys,
            args.k,
            args.min_novel_kmers,
            args.join_gap_kmers,
            args.flank,
        )
        intervals = [(s, e) for s, e in intervals if e - s >= args.min_segment_length]
        if not intervals:
            fully_redundant_records += 1
            rows.append((name, len(seq), informative, represented, fraction, "redundant", 0, 0))
            continue

        segmented_records += 1
        segmented_source_bases += len(seq)
        emitted = 0
        for idx, (start, end) in enumerate(intervals, 1):
            segment = seq[start:end]
            kept_segments.append((f"{name}.novel{idx} source={name}:{start}-{end}", segment))
            emitted += len(segment)
            novel_segment_bases += len(segment)
        rows.append((name, len(seq), informative, represented, fraction, "segmented", len(intervals), emitted))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        index = 0
        for source, records in (("backbone", backbone), ("recovery", kept_whole), ("novel", kept_segments)):
            for name, seq in records:
                index += 1
                write_record(out, f"{source}_{index:07d} original={name}", seq)

    with args.report.open("w") as out:
        out.write("contig\tlength\tinformative_kmers\trepresented_kmers\trepresented_fraction\taction\tsegments\temitted_bases\n")
        for row in rows:
            out.write("\t".join(map(str, row[:4])) + f"\t{row[4]:.6f}\t" + "\t".join(map(str, row[5:])) + "\n")

    stats = {
        "backbone_records": len(backbone),
        "backbone_bases": sum(len(seq) for _, seq in backbone),
        "recovery_records": len(recovery),
        "recovery_bases": sum(len(seq) for _, seq in recovery),
        "recovery_whole_records": len(kept_whole),
        "recovery_whole_bases": sum(len(seq) for _, seq in kept_whole),
        "recovery_segmented_records": segmented_records,
        "recovery_segmented_source_bases": segmented_source_bases,
        "novel_segment_records": len(kept_segments),
        "novel_segment_bases": novel_segment_bases,
        "fully_redundant_records": fully_redundant_records,
        "segment_fraction": args.segment_fraction,
        "k": args.k,
        "output_records": len(backbone) + len(kept_whole) + len(kept_segments),
        "output_bases": sum(len(seq) for _, seq in backbone) + sum(len(seq) for _, seq in kept_whole) + novel_segment_bases,
    }
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, sort_keys=True))


if __name__ == "__main__":
    main()
