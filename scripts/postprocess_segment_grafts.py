#!/usr/bin/env python3
"""Postprocess segment-level recovery with stitching before the 200 bp cutoff.

Short novel tails can be biologically useful even when the tail+anchor record is
<200 bp: if it has an exact unique overlap with a long backbone, it should be
grafted first and only then subjected to the final contig length filter. This
avoids padding each recovery segment with represented backbone merely to survive
early filtering.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(cmd: list[object]) -> None:
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--short-min-length", type=int, default=64)
    ap.add_argument("--final-min-length", type=int, default=200)
    ap.add_argument("--min-overlap", type=int, default=31)
    ap.add_argument("--overlap-margin", type=int, default=10)
    args = ap.parse_args()
    if args.short_min_length < args.min_overlap:
        raise SystemExit("short-min-length must be >= min-overlap")
    if args.final_min_length < args.short_min_length:
        raise SystemExit("final-min-length must be >= short-min-length")

    scripts = Path(__file__).resolve().parent
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    merged = out / "short_union.fasta"
    prefiltered = out / "short_noncontained.fasta"
    prestitched = out / "short_stitched.fasta"
    long_filtered = out / "long_noncontained.fasta"
    final = out / "primary_contigs.fasta"

    run([
        sys.executable,
        scripts / "merge_fasta_unique.py",
        merged,
        args.input,
        "--min-length",
        args.short_min_length,
    ])
    run([
        sys.executable,
        scripts / "filter_contained_fasta.py",
        merged,
        prefiltered,
        "--min-length",
        args.short_min_length,
        "--seed-k",
        21,
        "--window",
        12,
        "--candidate-minimizers",
        16,
        "--removed-tsv",
        out / "short_contained_removed.tsv",
        "--stats-json",
        out / "short_containment_stats.json",
    ])
    run([
        sys.executable,
        scripts / "stitch_exact_overlaps.py",
        prestitched,
        prefiltered,
        "--min-overlap",
        args.min_overlap,
        "--overlap-margin",
        args.overlap_margin,
        "--seed-length",
        args.min_overlap,
        "--max-seed-occurrences",
        64,
        "--min-length",
        args.short_min_length,
    ])
    run([
        sys.executable,
        scripts / "filter_contained_fasta.py",
        prestitched,
        long_filtered,
        "--min-length",
        args.final_min_length,
        "--seed-k",
        21,
        "--window",
        12,
        "--candidate-minimizers",
        16,
        "--removed-tsv",
        out / "long_contained_removed.tsv",
        "--stats-json",
        out / "long_containment_stats.json",
    ])
    run([
        sys.executable,
        scripts / "stitch_exact_overlaps.py",
        final,
        long_filtered,
        "--min-overlap",
        args.min_overlap,
        "--overlap-margin",
        args.overlap_margin,
        "--seed-length",
        args.min_overlap,
        "--max-seed-occurrences",
        64,
        "--min-length",
        args.final_min_length,
    ])

    def count(path: Path) -> tuple[int, int]:
        n = bases = 0
        current = 0
        for raw in path.read_text().splitlines():
            if raw.startswith(">"):
                if current:
                    n += 1
                    bases += current
                current = 0
            else:
                current += len(raw.strip())
        if current:
            n += 1
            bases += current
        return n, bases

    stats = {}
    for label, path in (
        ("input", args.input),
        ("short_prefiltered", prefiltered),
        ("short_stitched", prestitched),
        ("long_filtered", long_filtered),
        ("final", final),
    ):
        records, bases = count(path)
        stats[f"{label}_records"] = records
        stats[f"{label}_bases"] = bases
    (out / "segment_graft_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
