#!/usr/bin/env python3
"""Segment-level reference-free merge of a graph backbone and recovery contigs.

Recovery contigs that mostly duplicate the graph backbone are not discarded as
whole records. Instead, only novel k-mer intervals are retained, together with
the minimum represented flank needed for exact-overlap stitching. Terminal
novel sequence carries an anchor only on its represented side; internal novel
sequence carries one anchor on each represented side. This keeps the contiguity
benefit of recovery segments without re-emitting unnecessary backbone sequence.
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
    left_anchor: int,
    right_anchor: int,
    min_length: int,
) -> tuple[int, int]:
    """Expand only into represented flanks that are actually available.

    A terminal novel run has represented sequence on only one side, so adding
    an anchor on the terminal side can only duplicate sequence and cannot help
    an exact-overlap graft. Internal runs retain anchors on both sides.
    """
    start = max(0, start - left_anchor)
    end = min(seq_len, end + right_anchor)
    missing = min_length - (end - start)
    if missing <= 0:
        return start, end

    # If a very short novel segment still needs padding to survive the short
    # prefilter, consume represented flank only on sides that already have an
    # anchor. This avoids inventing a terminal duplicate merely to reach a size
    # threshold.
    while missing > 0 and (left_anchor > 0 or right_anchor > 0):
        changed = False
        if left_anchor > 0 and start > 0:
            take = min(start, max(1, (missing + 1) // 2))
            start -= take
            missing -= take
            changed = True
        if missing > 0 and right_anchor > 0 and end < seq_len:
            take = min(seq_len - end, missing)
            end += take
            missing -= take
            changed = True
        if not changed:
            break
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
    ap.add_argument("--stats-json", type=Path)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("--replace-fraction", type=float, default=0.85)
    ap.add_argument("--min-informative-kmers", type=int, default=20)
    ap.add_argument("--segment-anchor-bases", type=int, default=31)
    ap.add_argument("--min-novel-kmers", type=int, default=4)
    ap.add_argument("--merge-represented-gap-kmers", type=int, default=8)
    ap.add_argument("--max-novel-hole-kmers", type=int, default=2)
    ap.add_argument("--min-segment-length", type=int, default=31)
    args = ap.parse_args()

    if not 0.0 <= args.replace_fraction <= 1.0:
        raise SystemExit("replace-fraction must be in [0,1]")
    if args.k < 3:
        raise SystemExit("k must be >=3")
    if args.segment_anchor_bases < args.k:
        raise SystemExit("segment-anchor-bases must be >= k for exact-overlap grafting")
    if args.min_segment_length < args.k:
        raise SystemExit("min-segment-length must be >= k")

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
        "recovery_anchor_bases_retained": 0,
        "terminal_segment_records": 0,
        "internal_segment_records": 0,
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
                    0,
                    "full",
                )
            )
            continue

        fill_short_novel_holes(state, args.max_novel_hole_kmers)
        candidate_runs = novel_runs(state, args.merge_represented_gap_kmers)
        intervals: list[tuple[int, int, int, str]] = []
        state_len = len(state)
        for start_k, end_k in candidate_runs:
            novel_count = sum(
                1 for value in state[start_k:end_k] if value is not True
            )
            if novel_count < args.min_novel_kmers:
                continue
            base_start = start_k
            base_end = min(len(seq), end_k + args.k - 1)
            left_represented = start_k > 0 and state[start_k - 1] is True
            right_represented = end_k < state_len and state[end_k] is True
            left_anchor = args.segment_anchor_bases if left_represented else 0
            right_anchor = args.segment_anchor_bases if right_represented else 0
            expanded_start, expanded_end = expand_interval(
                base_start,
                base_end,
                len(seq),
                left_anchor,
                right_anchor,
                args.min_segment_length,
            )
            boundary = (
                "internal"
                if left_represented and right_represented
                else "terminal"
                if left_represented or right_represented
                else "unanchored"
            )
            intervals.append((expanded_start, expanded_end, novel_count, boundary))

        # Preserve boundary metadata while merging overlapping extraction windows.
        merged_coords = merge_intervals([(start, end) for start, end, _, _ in intervals])
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
                    0,
                    "none",
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
            represented_count = sum(
                1 for value in state[left_k:right_k] if value is True
            )
            # Classify the merged segment by whether it touches either recovery end.
            boundary = "terminal" if start == 0 or end == len(seq) else "internal"
            stats[f"{boundary}_segment_records"] += 1
            stats["recovery_novel_kmers_retained"] += novel_count
            stats["recovery_anchor_bases_retained"] += represented_count
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
                    represented_count,
                    boundary,
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
            "segment_length\tnovel_kmers\trepresented_segment_kmers\tboundary\n"
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
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, sort_keys=True))


if __name__ == "__main__":
    main()
