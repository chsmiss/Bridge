#!/usr/bin/env python3
"""Partition paired FASTQ reads into reference-genome bins using PAF evidence.

This is a diagnostic/oracle utility, not a production binning algorithm. Reference
identifiers are expected to be prefixed as ``BIN|original_contig``. Mate 1 and
mate 2 are aligned separately so evidence from reads sharing the same spot name
cannot be mixed. A pair is assigned to one bin when both mates support the same
unique bin, or when one mate is uniquely assigned and the other mate has no
conflicting confident assignment. Ambiguous/conflicting/unmapped pairs are kept
in a shared pool rather than discarded.
"""
from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


@dataclass(frozen=True)
class Hit:
    score: int
    mapq: int
    aligned_fraction: float


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode)


def normalize_name(text: str) -> str:
    name = text.strip()
    if name.startswith("@"):
        name = name[1:]
    name = name.split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def fastq_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with open_text(path, "r") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ: {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record in {path}")
            yield (
                header.rstrip("\n"),
                sequence.rstrip("\n"),
                plus.rstrip("\n"),
                quality.rstrip("\n"),
            )


def paf_hits(path: Path, delimiter: str) -> dict[str, dict[str, Hit]]:
    by_read: dict[str, dict[str, Hit]] = defaultdict(dict)
    with path.open() as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"invalid PAF line: {raw[:120]!r}")
            qname = normalize_name(fields[0])
            qlen = max(1, int(fields[1]))
            qstart = int(fields[2])
            qend = int(fields[3])
            genome = fields[5].split(delimiter, 1)[0]
            matches = int(fields[9])
            mapq = int(fields[11])
            score = matches
            for tag in fields[12:]:
                if tag.startswith("AS:i:"):
                    score = int(tag[5:])
                    break
            aligned_fraction = max(0, qend - qstart) / qlen
            hit = Hit(score=score, mapq=mapq, aligned_fraction=aligned_fraction)
            previous = by_read[qname].get(genome)
            if previous is None or (hit.score, hit.mapq, hit.aligned_fraction) > (
                previous.score,
                previous.mapq,
                previous.aligned_fraction,
            ):
                by_read[qname][genome] = hit
    return by_read


def unique_assignment(
    hits: dict[str, Hit] | None,
    min_mapq: int,
    min_fraction: float,
    score_margin: int,
) -> tuple[str | None, str]:
    if not hits:
        return None, "unmapped"
    candidates = [
        (genome, hit)
        for genome, hit in hits.items()
        if hit.mapq >= min_mapq and hit.aligned_fraction >= min_fraction
    ]
    if not candidates:
        return None, "weak"
    candidates.sort(
        key=lambda item: (item[1].score, item[1].mapq, item[1].aligned_fraction),
        reverse=True,
    )
    best_genome, best = candidates[0]
    if len(candidates) == 1:
        return best_genome, "unique"
    second = candidates[1][1]
    if best.score - second.score >= score_margin:
        return best_genome, "unique"
    return None, "ambiguous"


def safe_label(label: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("._")
    return cleaned or "unnamed"


def write_record(handle: TextIO, record: tuple[str, str, str, str]) -> None:
    header, sequence, plus, quality = record
    handle.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf1", required=True, type=Path)
    parser.add_argument("--paf2", required=True, type=Path)
    parser.add_argument("--read1", required=True, type=Path)
    parser.add_argument("--read2", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target-delimiter", default="|")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.70)
    parser.add_argument("--score-margin", type=int, default=10)
    parser.add_argument(
        "--strict-both-mates",
        action="store_true",
        help="Require both mates to be uniquely assigned to the same bin.",
    )
    args = parser.parse_args()

    if not 0 <= args.min_mapq <= 255:
        raise ValueError("min-mapq must be in 0..255")
    if not 0.0 <= args.min_aligned_fraction <= 1.0:
        raise ValueError("min-aligned-fraction must be in 0..1")
    if args.score_margin < 0:
        raise ValueError("score-margin must be non-negative")

    args.output.mkdir(parents=True, exist_ok=True)
    bins_dir = args.output / "bins"
    bins_dir.mkdir(exist_ok=True)
    left_hits = paf_hits(args.paf1, args.target_delimiter)
    right_hits = paf_hits(args.paf2, args.target_delimiter)

    handles: dict[str, tuple[TextIO, TextIO]] = {}
    counts: Counter[str] = Counter()
    reasons: Counter[str] = Counter()

    def get_handles(label: str) -> tuple[TextIO, TextIO]:
        if label not in handles:
            directory = args.output / "shared" if label == "shared" else bins_dir / safe_label(label)
            directory.mkdir(parents=True, exist_ok=True)
            handles[label] = (
                gzip.open(directory / "R1.fastq.gz", "wt"),
                gzip.open(directory / "R2.fastq.gz", "wt"),
            )
        return handles[label]

    assignment_path = args.output / "pair_assignments.tsv.gz"
    with gzip.open(assignment_path, "wt") as assignment_out:
        assignment_out.write("pair\tbin\tleft_status\tright_status\treason\n")
        left_iter = fastq_records(args.read1)
        right_iter = fastq_records(args.read2)
        pair_index = 0
        while True:
            try:
                left = next(left_iter)
            except StopIteration:
                left = None
            try:
                right = next(right_iter)
            except StopIteration:
                right = None
            if left is None and right is None:
                break
            if left is None or right is None:
                raise ValueError("paired FASTQ files contain different record counts")
            pair_index += 1
            left_name = normalize_name(left[0])
            right_name = normalize_name(right[0])
            if left_name != right_name:
                raise ValueError(
                    f"pair name mismatch at pair {pair_index}: {left_name!r} != {right_name!r}"
                )

            left_bin, left_status = unique_assignment(
                left_hits.get(left_name),
                args.min_mapq,
                args.min_aligned_fraction,
                args.score_margin,
            )
            right_bin, right_status = unique_assignment(
                right_hits.get(right_name),
                args.min_mapq,
                args.min_aligned_fraction,
                args.score_margin,
            )

            output_bin = "shared"
            reason = "ambiguous_or_unmapped"
            if left_bin is not None and right_bin is not None:
                if left_bin == right_bin:
                    output_bin = left_bin
                    reason = "both_mates_agree"
                else:
                    reason = "mate_conflict"
            elif not args.strict_both_mates:
                if left_bin is not None and right_bin is None and right_status in {"unmapped", "weak"}:
                    output_bin = left_bin
                    reason = "left_unique_mate_uninformative"
                elif right_bin is not None and left_bin is None and left_status in {"unmapped", "weak"}:
                    output_bin = right_bin
                    reason = "right_unique_mate_uninformative"

            out1, out2 = get_handles(output_bin)
            write_record(out1, left)
            write_record(out2, right)
            counts[output_bin] += 1
            reasons[reason] += 1
            assignment_out.write(
                f"{left_name}\t{output_bin}\t{left_status}\t{right_status}\t{reason}\n"
            )

    for out1, out2 in handles.values():
        out1.close()
        out2.close()

    manifest = args.output / "manifest.tsv"
    with manifest.open("w") as handle:
        handle.write("bin\tpairs\tread1\tread2\n")
        for label in sorted(counts, key=lambda value: (value == "shared", value)):
            directory = args.output / "shared" if label == "shared" else bins_dir / safe_label(label)
            handle.write(
                f"{label}\t{counts[label]}\t{directory / 'R1.fastq.gz'}\t{directory / 'R2.fastq.gz'}\n"
            )

    with (args.output / "assignment_summary.tsv").open("w") as handle:
        handle.write("category\tcount\n")
        for reason, count in sorted(reasons.items()):
            handle.write(f"{reason}\t{count}\n")
        handle.write(f"total_pairs\t{sum(counts.values())}\n")
        handle.write(f"reference_bins\t{sum(1 for key in counts if key != 'shared')}\n")
        handle.write(f"shared_pairs\t{counts['shared']}\n")

    print(manifest.read_text(), end="")
    print((args.output / "assignment_summary.tsv").read_text(), end="")


if __name__ == "__main__":
    main()
