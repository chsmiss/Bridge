#!/usr/bin/env python3
"""Stage30: challenge k31 path-cover continuity behind an immutable Stage24 seed.

The experiment deliberately separates candidate generation from production
acceptance.  k31 can use progressively more permissive path-cover settings,
but every final assembly is produced by the same Stage24 seed-lock500 rule.
Therefore an aggressive k31 candidate cannot replace or delete a validated
Stage24 contig; it can only contribute a unique reciprocal exact end extension.

This is reference-free.  Benchmark workflows use references only after the
assemblies have been generated.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import stage26_stage24_carry_forward as s26


CONFIGS = (
    {
        "name": "threaded_strict",
        "major_path_cover": False,
        "min_primary_support": 5,
        "primary_dominance": 0.75,
        "secondary_dominance": 0.25,
    },
    {
        "name": "major_strict",
        "major_path_cover": True,
        "min_primary_support": 5,
        "primary_dominance": 0.75,
        "secondary_dominance": 0.25,
    },
    {
        "name": "major_balanced",
        "major_path_cover": True,
        "min_primary_support": 4,
        "primary_dominance": 0.65,
        "secondary_dominance": 0.20,
    },
    {
        "name": "major_aggressive",
        "major_path_cover": True,
        "min_primary_support": 2,
        "primary_dominance": 0.55,
        "secondary_dominance": 0.10,
    },
)


def replace_option(command: list[str], flag: str, value: object) -> None:
    """Replace the value following one existing CLI flag in-place."""
    try:
        index = command.index(flag)
    except ValueError as exc:
        raise ValueError(f"missing expected assembler flag {flag}") from exc
    if index + 1 >= len(command):
        raise ValueError(f"assembler flag {flag} has no value")
    command[index + 1] = str(value)


def configured_command(
    *,
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    output: Path,
    threads: int,
    config: dict[str, object],
) -> list[str]:
    command = s26.seeded_assemble_command(
        bridgeasm,
        read1,
        read2,
        output,
        31,
        16,
        threads,
        major_path_cover=bool(config["major_path_cover"]),
    )
    replace_option(command, "--min-primary-support", config["min_primary_support"])
    replace_option(command, "--primary-dominance", config["primary_dominance"])
    replace_option(
        command,
        "--path-cover-secondary-dominance",
        config["secondary_dominance"],
    )
    return command


def validate_configs() -> None:
    names = [str(config["name"]) for config in CONFIGS]
    if len(names) != len(set(names)):
        raise ValueError("Stage30 config names must be unique")
    expected = {
        "threaded_strict": (False, 5, 0.75, 0.25),
        "major_strict": (True, 5, 0.75, 0.25),
        "major_balanced": (True, 4, 0.65, 0.20),
        "major_aggressive": (True, 2, 0.55, 0.10),
    }
    for config in CONFIGS:
        got = (
            bool(config["major_path_cover"]),
            int(config["min_primary_support"]),
            float(config["primary_dominance"]),
            float(config["secondary_dominance"]),
        )
        if got != expected[str(config["name"])]:
            raise ValueError(f"unexpected Stage30 config {config}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeasm", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--seed-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--prior-min-length", type=int, default=200)
    parser.add_argument("--seed-lock-min-overlap", type=int, default=500)
    parser.add_argument("--seed-lock-overlap-margin", type=int, default=30)
    parser.add_argument("--seed-lock-min-extension", type=int, default=20)
    parser.add_argument("--list-configs", action="store_true")
    args = parser.parse_args()

    validate_configs()
    if args.list_configs:
        print(json.dumps(CONFIGS, indent=2))
        return

    required = [args.bridgeasm, args.read1, args.read2, args.seed_fasta]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage30 inputs: " + ", ".join(missing))
    if args.prior_min_length < 1:
        raise SystemExit("--prior-min-length must be positive")
    if args.seed_lock_min_overlap < 31:
        raise SystemExit("--seed-lock-min-overlap must be >=31")

    args.output.mkdir(parents=True, exist_ok=True)
    virtual_r1 = args.output / "virtual_k31_R1.fastq.gz"
    virtual_r2 = args.output / "virtual_k31_R2.fastq.gz"
    virtual_records, virtual_bases = s26.carry.append_virtual_pairs(
        args.read1,
        args.read2,
        args.seed_fasta,
        virtual_r1,
        virtual_r2,
        args.prior_min_length,
    )

    scripts = Path(__file__).resolve().parent
    results: list[dict[str, object]] = []
    try:
        for config in CONFIGS:
            name = str(config["name"])
            root = args.output / name
            k31 = root / "k31"
            assemble_seconds = s26.carry.run(
                configured_command(
                    bridgeasm=args.bridgeasm,
                    read1=virtual_r1,
                    read2=virtual_r2,
                    output=k31,
                    threads=args.threads,
                    config=config,
                ),
                env=s26.clean_recovery_env(),
            )

            final_fasta = root / "primary_contigs.fasta"
            seed_lock_stats = root / "seed_lock_stats.json"
            seed_lock_seconds = s26.carry.run(
                s26.seed_lock_command(
                    scripts,
                    final_fasta,
                    args.seed_fasta,
                    k31 / "primary_contigs.fasta",
                    seed_lock_stats,
                    min_overlap=args.seed_lock_min_overlap,
                    overlap_margin=args.seed_lock_overlap_margin,
                    min_extension=args.seed_lock_min_extension,
                )
            )
            lock = json.loads(seed_lock_stats.read_text())
            results.append(
                {
                    **config,
                    "assemble_seconds": assemble_seconds,
                    "seed_lock_seconds": seed_lock_seconds,
                    "accepted_extensions": lock.get("accepted_extensions", 0),
                    "added_bp": lock.get("added_bp", 0),
                    "ambiguous_bridge_candidates_rejected": lock.get(
                        "ambiguous_bridge_candidates_rejected", 0
                    ),
                    "final_fasta": str(final_fasta),
                }
            )
    finally:
        virtual_r1.unlink(missing_ok=True)
        virtual_r2.unlink(missing_ok=True)

    manifest = {
        "pipeline": "stage30-seedlock-pathcover-sweep-v1",
        "seed_fasta": str(args.seed_fasta),
        "prior_min_length": args.prior_min_length,
        "virtual_records": virtual_records,
        "virtual_bases": virtual_bases,
        "seed_lock_min_overlap": args.seed_lock_min_overlap,
        "seed_lock_overlap_margin": args.seed_lock_overlap_margin,
        "seed_lock_min_extension": args.seed_lock_min_extension,
        "configs": results,
    }
    (args.output / "sweep_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
