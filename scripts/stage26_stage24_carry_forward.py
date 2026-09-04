#!/usr/bin/env python3
"""Stage26: preserve the Stage24 pair2 assembly while improving continuity.

Stage24 established a conservative singleton-rescue seed: target-k singleton
nodes are admitted only from previous-k topology with observed target-k edges
and support from at least two distinct physical fragments.  Stage26 treats that
assembly as immutable sequence evidence and carries it through fresh k31/k41/k55
assemblies using at most one virtual-fragment vote per prior contig.

The final catalog always includes the Stage24 seed before exact containment and
exact-overlap stitching.  Therefore continuity work is not allowed to improve
N50 by simply deleting Stage24 sequence.  The workflow separately measures
Stage24-seed k21 retention to make this invariant visible.

This pipeline is reference-free.  References are used only by benchmark
workflows to evaluate GF and assembly correctness.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import bridge_carry_forward_pipeline as carry


def seeded_assemble_command(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    output: Path,
    k: int,
    mercy: int,
    threads: int,
    *,
    major_path_cover: bool,
) -> list[str]:
    """Build the evidence-aware assembler command for one carry-forward stage."""
    command = carry.assemble_command(
        bridgeasm,
        read1,
        read2,
        output,
        k,
        mercy,
        threads,
    )
    if not major_path_cover:
        command.remove("--major-path-cover")
    return command


def clean_recovery_env() -> dict[str, str]:
    """Keep Stage26 from silently re-enabling older singleton experiments."""
    env = dict(os.environ)
    for key in (
        "BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION",
        "BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY",
        "BRIDGEASM_MATE_TERMINAL_MERCY_KMERS",
    ):
        env.pop(key, None)
    return env


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeasm", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--seed-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--prior-min-length", type=int, default=200)
    parser.add_argument("--major-path-cover", action="store_true")
    parser.add_argument("--stitch-min-overlap", type=int, default=40)
    parser.add_argument("--stitch-overlap-margin", type=int, default=10)
    args = parser.parse_args()

    required = [args.bridgeasm, args.read1, args.read2, args.seed_fasta]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage26 inputs: " + ", ".join(missing))
    if args.prior_min_length < 1:
        raise SystemExit("--prior-min-length must be positive")

    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    timings: dict[str, float] = {}
    virtual: dict[str, dict[str, int]] = {}
    stage_dirs: list[Path] = []

    previous_fasta = args.seed_fasta
    stage_specs = [(31, 16), (41, 12), (55, 8)]
    for k, mercy in stage_specs:
        augmented_r1 = output / f"virtual_k{k}_R1.fastq.gz"
        augmented_r2 = output / f"virtual_k{k}_R2.fastq.gz"
        count, bases = carry.append_virtual_pairs(
            args.read1,
            args.read2,
            previous_fasta,
            augmented_r1,
            augmented_r2,
            args.prior_min_length,
        )
        virtual[f"k{k}"] = {"records": count, "bases": bases}

        stage = output / f"k{k}"
        timings[f"k{k}"] = carry.run(
            seeded_assemble_command(
                args.bridgeasm,
                augmented_r1,
                augmented_r2,
                stage,
                k,
                mercy,
                args.threads,
                major_path_cover=args.major_path_cover,
            ),
            env=clean_recovery_env(),
        )
        augmented_r1.unlink(missing_ok=True)
        augmented_r2.unlink(missing_ok=True)
        previous_fasta = stage / "primary_contigs.fasta"
        stage_dirs.append(stage)

    # The seed is deliberately included in the final union.  Higher-k stages
    # are allowed to supersede/contain it, but never to erase it merely because
    # a rescued low-depth path cannot survive at a larger k.
    candidates: list[Path] = [args.seed_fasta]
    for stage in stage_dirs:
        candidates.append(stage / "primary_contigs.fasta")
        candidates.append(stage / "haplotigs.fasta")

    union = output / "cross_k_exact_union.fasta"
    timings["union"] = carry.run(
        [
            sys.executable,
            str(scripts / "merge_fasta_unique.py"),
            str(union),
            *map(str, candidates),
            "--min-length",
            "200",
        ]
    )

    filtered = output / "cross_k_noncontained.fasta"
    timings["containment"] = carry.run(
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
    timings["stitch"] = carry.run(
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
        "pipeline": "stage26-stage24-pair2-carry-forward-v1",
        "seed_fasta": str(args.seed_fasta),
        "prior_min_length": args.prior_min_length,
        "major_path_cover": args.major_path_cover,
        "virtual_reads": virtual,
        "stages": [31, 41, 55],
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
