#!/usr/bin/env python3
"""Segment-level reference-free merge of a graph backbone and recovery contigs.

Recovery contigs that mostly duplicate the graph backbone are not discarded as
whole records. Instead, only novel k-mer intervals are retained, together with
short exact-sequence anchors from the represented flanks. This preserves novel
extensions/bridges for downstream exact-overlap stitching while avoiding the
large duplication penalty of emitting the complete overlapping recovery contig.
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


def canonical_kmers(seq: str, k: int):
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


def fill_short_novel_holes(state: list[bool | None], max_hole: int) -> None:
    """Suppress tiny novel islands surrounded by represented k-mers."""
    if max_hole <= 0:
        return
    i = 0
    while i < len(state):
        if state[i] is True:
            i += 1
            continue
        start = i
        while i < len(state) and state[i] is not True:
            i += 1
        end = i
        if (
            start > 0
            and end < len(state)
            and state[start - 1] is True
            and state[end] is True
            and end - start <= max_hole
        ):
            for j in range(start, end):
                state[j] = True


def novel_runs(state: list[bool | None], merge_gap: int) -> list[tuple[int, int]]:
    """Return half-open k-mer-coordinate runs not represented by backbone."""
    runs: list[tuple[int, int]] = []
    i = 0
    while i < len(state):
        if state[i] is True:
            i += 1
            continue
        start = i
        while i < len(state) and state[i] is not True:
            i += 1
        runs.append((start, i))
    if not runs or merge_gap <= 0:
        return runs
    merged = [runs[0]]
    for start, end in runs[1:]:
        prev_start, prev_end = merged[-1]
        if start - prev_end <= merge_gap:
            merged[-1] = (prev_start, end)
        else:
            merged.append((start, end))
    return merged


def expand_interval(
    start: int,
    end: int,
    seq_len: int,
    anchor_bases: int,
    min_length: int,
) -> tuple[int, int]:
    start = max(0, start - anchor_bases)
    end = min(seq_len, end + anchor_bases)
    missing = min_length - (end - start)
    if missing > 0:
        left_room = start
        right_room = seq_len - end
        take_left = min(left_room, (missing + 1) // 2)
        start -= take_left
        missing -= take_left
        take_right = min(right_room, missing)
        end += take_right
        missing -= take_right
        if missing > 0:
            take_left = min(start, missing)
            start -= take_left
    return start, end


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", type=Path, required=True)
    ap.add_argument("--recovery", type=Path, required=True)
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("--replace-fraction", type=float, default=0.85)
    ap.add_argument("--min-informative-kmers", type=int, default=20)
    ap.add_argument("--segment-anchor-bases", type=int, default=96)
    ap.add_argument("--min-novel-kmers", type=int, default=4)
    ap.add_argument("--merge-represented-gap-kmers", type=int, default=8)
    ap.add_argument("--max-novel-hole-kmers", type=int, default=2)
    ap.add_argument("--min-segment-length", type=int, default=200)
    args = ap.parse_args()

    if not 0.0 <= args.replace_fraction <= 1.0:
        raise SystemExit("replace-fraction must be in [0,1]")
    if args.k < 3:
        raise SystemExit("k must be >=3")

    backbone = list(fasta(args.backbone))
    recovery = list(fasta(args.recovery))
    backbone_keys: set[int] = set()
    for _, seq in backbone:
        backbone_keys.update(key for _, key in canonical_kmers(seq, args.k))

    emitted: list[tuple[str, str, str]] = []
    rows: list[tuple[object, ...]] = []
    stats = {
        "backbone_records": len(backbone),
        "backbone_bases": sum(len(seq) for _, seq in backbone),
        "recovery_records": len(recovery),
        "recovery_bases": sum(len(seq) for _, seq in recovery),
        "recovery_kept_full_records": 0,
        "recovery_kept_full_bases": 0,
        "recovery_segmented_records": 0,
        "recovery_segment_records": 0,
        "recovery_segment_bases": 0,
        "recovery_fully_replaced_records": 0,
        "recovery_fully_replaced_bases": 0,
        "recovery_novel_kmers_retained": 0,
        "replace_fraction": args.replace_fraction,
        "k": args.k,
        "segment_anchor_bases": args.segment_anchor_bases,
    }

    for name, seq in recovery:
        state: list[bool | None] = [None] * max(0, len(seq) - args.k + 1)
        informative = represented = 0
        for pos, key in canonical_kmers(seq, args.k):
            hit = key in backbone_keys
            state[pos] = hit
            informative += 1
            represented += int(hit)
        fraction = represented / max(1, informative)

        if informative < args.min_informative_kmers or fraction < args.replace_fraction:
            emitted.append(("recovery_full", name, seq))
            stats["recovery_kept_full_records"] += 1
            stats["recovery_kept_full_bases"] += len(seq)
            rows.append(
                (
                    name,
                    len(seq),
                    informative,
                    represented,
                    fraction,
                    "keep_full",
                    0,
                    len(seq),
                    len(seq),
                    0,
                )
            )
            continue

        fill_short_novel_holes(state, args.max_novel_hole_kmers)
        candidate_runs = novel_runs(state, args.merge_represented_gap_kmers)
        intervals: list[tuple[int, int]] = []
        for start_k, end_k in candidate_runs:
            novel_count = sum(
                1 for value in state[start_k:end_k] if value is not True
            )
            if novel_count < args.min_novel_kmers:
                continue
            base_start = start_k
            base_end = min(len(seq), end_k + args.k - 1)
            intervals.append(
                expand_interval(
                    base_start,
                    base_end,
                    len(seq),
                    args.segment_anchor_bases,
                    args.min_segment_length,
                )
            )

        merged_coords = merge_intervals(intervals)
        if not merged_coords:
            stats["recovery_fully_replaced_records"] += 1
            stats["recovery_fully_replaced_bases"] += len(seq)
            rows.append(
                (
                    name,
                    len(seq),
                    informative,
                    represented,
                    fraction,
                    "replace_full",
                    0,
                    0,
                    0,
                    0,
                )
            )
            continue

        stats["recovery_segmented_records"] += 1
        for seg_idx, (start, end) in enumerate(merged_coords, 1):
            segment = seq[start:end]
            left_k = max(0, start)
            right_k = min(len(state), max(left_k, end - args.k + 1))
            novel_count = sum(
                1 for value in state[left_k:right_k] if value is not True
            )
            stats["recovery_novel_kmers_retained"] += novel_count
            stats["recovery_segment_records"] += 1
            stats["recovery_segment_bases"] += len(segment)
            emitted.append((f"recovery_segment_{seg_idx}", name, segment))
            rows.append(
                (
                    name,
                    len(seq),
                    informative,
                    represented,
                    fraction,
                    "keep_segment",
                    start,
                    end,
                    len(segment),
                    novel_count,
                )
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        index = 0
        for name, seq in backbone:
            index += 1
            out.write(f">backbone_{index:07d} original={name} len={len(seq)}\n")
            for start in range(0, len(seq), 80):
                out.write(seq[start : start + 80] + "\n")
        for source, name, seq in emitted:
            index += 1
            out.write(f">{source}_{index:07d} original={name} len={len(seq)}\n")
            for start in range(0, len(seq), 80):
                out.write(seq[start : start + 80] + "\n")

    with args.report.open("w") as out:
        out.write(
            "contig\tlength\tinformative_kmers\trepresented_kmers\t"
            "represented_fraction\taction\tsegment_start\tsegment_end\t"
            "segment_length\tnovel_kmers\n"
        )
        for row in rows:
            out.write(
                "\t".join(map(str, row[:4]))
                + f"\t{row[4]:.6f}\t"
                + "\t".join(map(str, row[5:]))
                + "\n"
            )

    stats["output_records"] = len(backbone) + len(emitted)
    stats["output_bases"] = stats["backbone_bases"] + sum(
        len(seq) for _, _, seq in emitted
    )
    print(json.dumps(stats, sort_keys=True))


if __name__ == "__main__":
    main()
