#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import stage26_stage24_carry_forward as s26


def test_major_path_cover_toggle() -> None:
    common = dict(
        bridgeasm=Path("bridgeasm"),
        read1=Path("r1.fq.gz"),
        read2=Path("r2.fq.gz"),
        output=Path("out"),
        k=31,
        mercy=16,
        threads=2,
    )
    conservative = s26.seeded_assemble_command(**common, major_path_cover=False)
    major = s26.seeded_assemble_command(**common, major_path_cover=True)
    assert "--threaded-path-cover" in conservative
    assert "--major-path-cover" not in conservative
    assert "--major-path-cover" in major
    assert conservative.count("--major-path-cover") == 0
    assert major.count("--major-path-cover") == 1


def test_old_singleton_experiments_are_removed_from_env() -> None:
    keys = (
        "BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION",
        "BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY",
        "BRIDGEASM_MATE_TERMINAL_MERCY_KMERS",
    )
    old = {key: os.environ.get(key) for key in keys}
    try:
        for key in keys:
            os.environ[key] = "sentinel"
        env = s26.clean_recovery_env()
        for key in keys:
            assert key not in env
    finally:
        for key, value in old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def test_seed_lock_is_the_default_stage_plan() -> None:
    assert s26.carry_stage_specs(legacy_multik=False) == [(31, 16)]
    assert s26.carry_stage_specs(legacy_multik=True) == [
        (31, 16),
        (41, 12),
        (55, 8),
    ]
    assert s26.DEFAULT_SEED_LOCK_MIN_OVERLAP == 500
    assert s26.DEFAULT_SEED_LOCK_OVERLAP_MARGIN == 30
    assert s26.DEFAULT_SEED_LOCK_MIN_EXTENSION == 20


def test_seed_lock_command_uses_k31_primary_and_immutable_seed() -> None:
    command = s26.seed_lock_command(
        Path("scripts"),
        Path("out/primary_contigs.fasta"),
        Path("seed.fasta"),
        Path("out/k31/primary_contigs.fasta"),
        Path("out/seed_lock_stats.json"),
        min_overlap=500,
        overlap_margin=30,
        min_extension=20,
    )
    joined = " ".join(map(str, command))
    assert "seed_locked_extensions.py" in joined
    assert "seed.fasta" in command
    assert "out/k31/primary_contigs.fasta" in command
    assert "--min-overlap" in command
    assert command[command.index("--min-overlap") + 1] == "500"
    assert command[command.index("--overlap-margin") + 1] == "30"
    assert command[command.index("--min-extension") + 1] == "20"
    assert "stitch_exact_overlaps.py" not in joined
    assert "filter_contained_fasta.py" not in joined


def main() -> None:
    test_major_path_cover_toggle()
    test_old_singleton_experiments_are_removed_from_env()
    test_seed_lock_is_the_default_stage_plan()
    test_seed_lock_command_uses_k31_primary_and_immutable_seed()
    print("stage26 carry-forward tests: passed")


if __name__ == "__main__":
    main()
