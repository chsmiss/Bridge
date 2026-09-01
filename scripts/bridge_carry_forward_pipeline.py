#!/usr/bin/env python3
"""Iterative multi-k BridgeAsm prototype with previous-k virtual long reads.

This is deliberately reference-free. The previous stage's sufficiently long
primary contigs are appended to read1 as high-quality synthetic reads. Matching
read2 records contain only N and therefore contribute no k-mers/threading and
cannot create a false mate bridge. One synthetic fragment contributes only one
unit of fragment support, so a prior path still needs support from the original
reads to cross the production min-count/min-fragment gates.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Iterator


def fasta_records(path: Path) -> Iterator[tuple[str, bytes]]:
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


def append_virtual_pairs(
    source_r1: Path,
    source_r2: Path,
    prior_fasta: Path,
    output_r1: Path,
    output_r2: Path,
    min_length: int,
) -> tuple[int, int]:
    output_r1.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_r1, output_r1)
    shutil.copyfile(source_r2, output_r2)
    kept = 0
    bases = 0
    # Appending creates an additional gzip member; flate2::MultiGzDecoder reads
    # concatenated members transparently.
    with gzip.open(output_r1, "ab", compresslevel=1) as left, gzip.open(
        output_r2, "ab", compresslevel=1
    ) as right:
        for index, (_header, sequence) in enumerate(fasta_records(prior_fasta), 1):
            if len(sequence) < min_length or set(sequence) - set(b"ACGT"):
                continue
            kept += 1
            bases += len(sequence)
            left.write(f"@bridge_prior_{index}\n".encode())
            left.write(sequence + b"\n+\n")
            left.write(b"I" * len(sequence) + b"\n")
            right.write(f"@bridge_prior_{index}\nN\n+\n!\n".encode())
    return kept, bases


def run(command: list[str], *, env: dict[str, str] | None = None) -> float:
    print("+", " ".join(command), flush=True)
    started = time.monotonic()
    subprocess.run(command, check=True, env=env)
    return time.monotonic() - started


def assemble_command(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    output: Path,
    k: int,
    mercy: int,
    threads: int,
) -> list[str]:
    return [
        str(bridgeasm),
        "assemble",
        "-1",
        str(read1),
        "-2",
        str(read2),
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
        "200",
        "--threads",
        str(threads),
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeasm", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--prior-min-length", type=int, default=500)
    parser.add_argument("--singleton-fraction", type=float, default=0.60)
    parser.add_argument("--mate-terminal-mercy", type=int, default=96)
    parser.add_argument("--stitch-min-overlap", type=int, default=40)
    parser.add_argument("--stitch-overlap-margin", type=int, default=10)
    args = parser.parse_args()

    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    timings: dict[str, float] = {}
    virtual: dict[str, dict[str, int]] = {}

    k21 = output / "k21"
    env = dict(**__import__("os").environ)
    env["BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION"] = str(args.singleton_fraction)
    env["BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY"] = "30"
    env["BRIDGEASM_MATE_TERMINAL_MERCY_KMERS"] = str(args.mate_terminal_mercy)
    timings["k21"] = run(
        assemble_command(args.bridgeasm, args.read1, args.read2, k21, 21, 24, args.threads),
        env=env,
    )

    previous = k21
    stage_specs = [(31, 16), (41, 12), (55, 8)]
    stage_dirs = [k21]
    for k, mercy in stage_specs:
        augmented_r1 = output / f"virtual_k{k}_R1.fastq.gz"
        augmented_r2 = output / f"virtual_k{k}_R2.fastq.gz"
        count, bases = append_virtual_pairs(
            args.read1,
            args.read2,
            previous / "primary_contigs.fasta",
            augmented_r1,
            augmented_r2,
            args.prior_min_length,
        )
        virtual[f"k{k}"] = {"records": count, "bases": bases}
        stage = output / f"k{k}"
        clean_env = dict(**__import__("os").environ)
        clean_env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION", None)
        clean_env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY", None)
        clean_env.pop("BRIDGEASM_MATE_TERMINAL_MERCY_KMERS", None)
        timings[f"k{k}"] = run(
            assemble_command(
                args.bridgeasm,
                augmented_r1,
                augmented_r2,
                stage,
                k,
                mercy,
                args.threads,
            ),
            env=clean_env,
        )
        augmented_r1.unlink(missing_ok=True)
        augmented_r2.unlink(missing_ok=True)
        previous = stage
        stage_dirs.append(stage)

    candidates: list[Path] = []
    for stage in stage_dirs:
        candidates.append(stage / "primary_contigs.fasta")
        candidates.append(stage / "haplotigs.fasta")
    union = output / "cross_k_exact_union.fasta"
    timings["union"] = run(
        [sys.executable, str(scripts / "merge_fasta_unique.py"), str(union), *map(str, candidates), "--min-length", "200"]
    )
    filtered = output / "cross_k_noncontained.fasta"
    timings["containment"] = run(
        [
            sys.executable,
            str(scripts / "filter_contained_fasta.py"),
            str(union),
            str(filtered),
            "--min-length",
            "200",
            "--seed-k",
            "21",
            "--window",
            "12",
            "--candidate-minimizers",
            "16",
            "--stats-json",
            str(output / "containment_stats.json"),
        ]
    )
    final_fasta = output / "primary_contigs.fasta"
    timings["stitch"] = run(
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
            "200",
        ]
    )

    manifest = {
        "pipeline": "bridge-carry-forward-v1",
        "prior_min_length": args.prior_min_length,
        "virtual_reads": virtual,
        "stages": [21, 31, 41, 55],
        "stitch_min_overlap": args.stitch_min_overlap,
        "stitch_overlap_margin": args.stitch_overlap_margin,
        "timings_seconds": timings,
        "total_seconds": sum(timings.values()),
        "final_fasta": str(final_fasta),
    }
    (output / "pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
