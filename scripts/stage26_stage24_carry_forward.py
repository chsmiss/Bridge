#!/usr/bin/env python3
"""Stage26/29: preserve the Stage24 pair2 assembly while improving continuity.

The production default is deliberately conservative:

* Stage24 pair2 is immutable sequence evidence.
* It is carried into one fresh k31 assembly using at most one virtual-fragment
  vote per prior contig.
* Only k31 primary contigs are considered for continuity recovery.
* A higher-k contig may extend one Stage24 contig end only through a long exact
  overlap that is the unique reciprocal best match.
* A candidate matching multiple physical Stage24 ends is rejected rather than
  used as a bridge; such bridges require independent read/mate evidence.

Stage28 showed that a 500-bp seed-lock overlap preserves every Stage24 canonical
k21 while improving N50/GF with near-baseline misassembly counts.  The older
k31->k41->k55 union/containment/stitch pipeline remains available only through
``--legacy-multik`` for reproducibility and diagnostics.

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


DEFAULT_SEED_LOCK_MIN_OVERLAP = 500
DEFAULT_SEED_LOCK_OVERLAP_MARGIN = 30
DEFAULT_SEED_LOCK_MIN_EXTENSION = 20


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


def carry_stage_specs(*, legacy_multik: bool) -> list[tuple[int, int]]:
    """Return k/mercy stages for the selected carry-forward mode."""
    if legacy_multik:
        return [(31, 16), (41, 12), (55, 8)]
    return [(31, 16)]


def seed_lock_command(
    scripts: Path,
    final_fasta: Path,
    seed_fasta: Path,
    k31_primary: Path,
    stats_json: Path,
    *,
    min_overlap: int,
    overlap_margin: int,
    min_extension: int,
) -> list[str]:
    """Build the immutable-seed, one-ended extension command."""
    return [
        sys.executable,
        str(scripts / "seed_locked_extensions.py"),
        str(final_fasta),
        str(seed_fasta),
        str(k31_primary),
        "--min-overlap",
        str(min_overlap),
        "--overlap-margin",
        str(overlap_margin),
        "--seed-length",
        "31",
        "--min-extension",
        str(min_extension),
        "--max-seed-occurrences",
        "64",
        "--stats-json",
        str(stats_json),
    ]


def run_legacy_multik(
    *,
    scripts: Path,
    output: Path,
    seed_fasta: Path,
    stage_dirs: list[Path],
    timings: dict[str, float],
    stitch_min_overlap: int,
    stitch_overlap_margin: int,
) -> Path:
    """Reproduce the pre-Stage29 cross-k union/containment/stitch behavior."""
    candidates: list[Path] = [seed_fasta]
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
            str(stitch_min_overlap),
            "--overlap-margin",
            str(stitch_overlap_margin),
            "--seed-length",
            "31",
            "--max-seed-occurrences",
            "64",
            "--min-length",
            "200",
        ]
    )
    return final_fasta


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
    parser.add_argument(
        "--legacy-multik",
        action="store_true",
        help="reproduce the old k31/k41/k55 union/containment/stitch pipeline",
    )
    parser.add_argument(
        "--seed-lock-min-overlap",
        type=int,
        default=DEFAULT_SEED_LOCK_MIN_OVERLAP,
    )
    parser.add_argument(
        "--seed-lock-overlap-margin",
        type=int,
        default=DEFAULT_SEED_LOCK_OVERLAP_MARGIN,
    )
    parser.add_argument(
        "--seed-lock-min-extension",
        type=int,
        default=DEFAULT_SEED_LOCK_MIN_EXTENSION,
    )
    parser.add_argument("--stitch-min-overlap", type=int, default=40)
    parser.add_argument("--stitch-overlap-margin", type=int, default=10)
    args = parser.parse_args()

    required = [args.bridgeasm, args.read1, args.read2, args.seed_fasta]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage26 inputs: " + ", ".join(missing))
    if args.prior_min_length < 1:
        raise SystemExit("--prior-min-length must be positive")
    if args.seed_lock_min_overlap < 31:
        raise SystemExit("--seed-lock-min-overlap must be >= 31")
    if args.seed_lock_overlap_margin < 0:
        raise SystemExit("--seed-lock-overlap-margin must be non-negative")
    if args.seed_lock_min_extension < 1:
        raise SystemExit("--seed-lock-min-extension must be positive")

    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    timings: dict[str, float] = {}
    virtual: dict[str, dict[str, int]] = {}
    stage_dirs: list[Path] = []

    previous_fasta = args.seed_fasta
    stage_specs = carry_stage_specs(legacy_multik=args.legacy_multik)
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

    if args.legacy_multik:
        final_fasta = run_legacy_multik(
            scripts=scripts,
            output=output,
            seed_fasta=args.seed_fasta,
            stage_dirs=stage_dirs,
            timings=timings,
            stitch_min_overlap=args.stitch_min_overlap,
            stitch_overlap_margin=args.stitch_overlap_margin,
        )
        mode = "legacy-multik"
        seed_lock_stats: dict[str, object] | None = None
    else:
        final_fasta = output / "primary_contigs.fasta"
        seed_lock_stats_path = output / "seed_lock_stats.json"
        timings["seed_lock"] = carry.run(
            seed_lock_command(
                scripts,
                final_fasta,
                args.seed_fasta,
                stage_dirs[0] / "primary_contigs.fasta",
                seed_lock_stats_path,
                min_overlap=args.seed_lock_min_overlap,
                overlap_margin=args.seed_lock_overlap_margin,
                min_extension=args.seed_lock_min_extension,
            )
        )
        seed_lock_stats = json.loads(seed_lock_stats_path.read_text())
        mode = "seed-lock-k31"

    manifest = {
        "pipeline": "stage26-stage24-pair2-carry-forward-v2",
        "mode": mode,
        "seed_fasta": str(args.seed_fasta),
        "prior_min_length": args.prior_min_length,
        "major_path_cover": args.major_path_cover,
        "virtual_reads": virtual,
        "stages": [k for k, _mercy in stage_specs],
        "seed_lock_min_overlap": None
        if args.legacy_multik
        else args.seed_lock_min_overlap,
        "seed_lock_overlap_margin": None
        if args.legacy_multik
        else args.seed_lock_overlap_margin,
        "seed_lock_min_extension": None
        if args.legacy_multik
        else args.seed_lock_min_extension,
        "seed_lock_stats": seed_lock_stats,
        "stitch_min_overlap": args.stitch_min_overlap if args.legacy_multik else None,
        "stitch_overlap_margin": args.stitch_overlap_margin
        if args.legacy_multik
        else None,
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
