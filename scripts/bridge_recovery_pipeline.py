#!/usr/bin/env python3
"""Reference-free BridgeAsm recovery/continuity pipeline.

The pipeline deliberately separates responsibilities across k values:

* k=21 keeps low-depth sequence with singleton-island plus mate-terminal rescue.
* k=31, k=41 and k=55 provide progressively more specific paths.
* all four assemblies enable same-read triplet and major-path evidence.
* physically-flanked, read-supported haplotigs are preserved as strain paths.
* exact reverse-complement-aware deduplication/containment filtering preserves
  complementary sequence without duplicating records.
* reciprocal unique exact suffix/prefix overlaps carry compatible paths across
  k values into longer final contigs.

No reference sequence, taxonomic label, or benchmark truth is used by any
assembly or merge decision.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time


def run(command: list[str], *, env: dict[str, str] | None = None) -> float:
    print("+", " ".join(command), flush=True)
    started = time.monotonic()
    subprocess.run(command, check=True, env=env)
    return time.monotonic() - started


def bridge_command(
    bridgeasm: Path,
    read1: Path,
    read2: Path | None,
    output: Path,
    k: int,
    mercy: int,
    threads: int,
    min_contig_length: int,
) -> list[str]:
    command = [
        str(bridgeasm),
        "assemble",
        "-1",
        str(read1),
        "-o",
        str(output),
        "-k",
        str(k),
        "--min-count",
        "2",
        "--mercy-max-kmers",
        str(mercy),
        "--mercy-min-support",
        "1",
        "--mercy-min-quality",
        "25",
        "--min-read-support",
        "2",
        "--min-pair-support",
        "2",
        "--min-primary-support",
        "5",
        "--primary-dominance",
        "0.75",
        "--threaded-path-cover",
        "--major-path-cover",
        "--path-cover-secondary-dominance",
        "0.25",
        "--min-contig-length",
        str(min_contig_length),
        "--threads",
        str(threads),
    ]
    if read2 is not None:
        command[4:4] = ["-2", str(read2)]
    return command


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeasm", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--min-contig-length", type=int, default=200)
    parser.add_argument("--singleton-fraction", type=float, default=0.60)
    parser.add_argument("--singleton-quality", type=float, default=30.0)
    parser.add_argument("--mate-terminal-mercy", type=int, default=96)
    parser.add_argument("--stitch-min-overlap", type=int, default=80)
    parser.add_argument("--stitch-overlap-margin", type=int, default=20)
    args = parser.parse_args()

    if args.threads <= 0:
        raise SystemExit("threads must be positive")
    if not 0.0 <= args.singleton_fraction <= 1.0:
        raise SystemExit("singleton fraction must be in [0,1]")
    if args.mate_terminal_mercy < 0:
        raise SystemExit("mate terminal mercy must be non-negative")
    if args.stitch_min_overlap < 31:
        raise SystemExit("stitch minimum overlap must be >=31")

    pipeline_started = time.monotonic()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    timings: dict[str, float] = {}

    stages = [
        ("k21_recall", 21, 24),
        ("k31_resolve", 31, 16),
        ("k41_resolve", 41, 12),
        ("k55_resolve", 55, 8),
    ]
    candidates: list[Path] = []
    for name, k, mercy in stages:
        stage_dir = output / name
        if stage_dir.exists():
            shutil.rmtree(stage_dir)
        env = os.environ.copy()
        if k == 21:
            env["BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION"] = str(
                args.singleton_fraction
            )
            env["BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY"] = str(
                args.singleton_quality
            )
            env["BRIDGEASM_MATE_TERMINAL_MERCY_KMERS"] = str(
                args.mate_terminal_mercy
            )
        else:
            env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION", None)
            env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY", None)
            env.pop("BRIDGEASM_MATE_TERMINAL_MERCY_KMERS", None)
        timings[name] = run(
            bridge_command(
                args.bridgeasm,
                args.read1,
                args.read2,
                stage_dir,
                k,
                mercy,
                args.threads,
                args.min_contig_length,
            ),
            env=env,
        )
        candidates.append(stage_dir / "primary_contigs.fasta")
        # Haplotigs are emitted only for bubbles that have unique flanks and
        # sufficient direct read support on both sides. Keeping them here is a
        # conservative way to preserve strain-specific paths that the single
        # primary path necessarily omits.
        candidates.append(stage_dir / "haplotigs.fasta")

    merged = output / "cross_k_exact_union.fasta"
    timings["exact_union"] = run(
        [
            sys.executable,
            str(scripts / "merge_fasta_unique.py"),
            str(merged),
            *map(str, candidates),
            "--min-length",
            str(args.min_contig_length),
        ]
    )

    filtered = output / "cross_k_noncontained.fasta"
    timings["containment_filter"] = run(
        [
            sys.executable,
            str(scripts / "filter_contained_fasta.py"),
            str(merged),
            str(filtered),
            "--min-length",
            str(args.min_contig_length),
            "--seed-k",
            "21",
            "--window",
            "12",
            "--candidate-minimizers",
            "16",
            "--removed-tsv",
            str(output / "contained_removed.tsv"),
            "--stats-json",
            str(output / "containment_stats.json"),
        ]
    )

    final_fasta = output / "primary_contigs.fasta"
    timings["exact_stitch"] = run(
        [
            sys.executable,
            str(scripts / "stitch_exact_overlaps.py"),
            str(final_fasta),
            str(filtered),
            "--min-overlap",
            str(args.stitch_min_overlap),
            "--overlap-margin",
            str(args.stitch_overlap_margin),
            "--seed-length",
            "31",
            "--max-seed-occurrences",
            "64",
            "--min-length",
            str(args.min_contig_length),
        ]
    )

    child_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest = {
        "pipeline": "bridge-recovery-v2",
        "read1": str(args.read1),
        "read2": str(args.read2) if args.read2 else None,
        "stages": [
            {"name": name, "k": k, "mercy_max_kmers": mercy}
            for name, k, mercy in stages
        ],
        "singleton_fraction_k21": args.singleton_fraction,
        "singleton_quality_k21": args.singleton_quality,
        "mate_terminal_mercy_k21": args.mate_terminal_mercy,
        "preserve_physically_flanked_haplotigs": True,
        "threaded_path_cover": True,
        "major_path_cover": True,
        "path_cover_secondary_dominance": 0.25,
        "stitch_min_overlap": args.stitch_min_overlap,
        "stitch_overlap_margin": args.stitch_overlap_margin,
        "final_fasta": str(final_fasta),
        "timings_seconds": timings,
        "total_stage_seconds": sum(timings.values()),
        "wall_seconds": time.monotonic() - pipeline_started,
        "peak_child_rss_mib": child_usage.ru_maxrss / 1024.0,
    }
    (output / "pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
