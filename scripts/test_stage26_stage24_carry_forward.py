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


def main() -> None:
    test_major_path_cover_toggle()
    test_old_singleton_experiments_are_removed_from_env()
    print("stage26 carry-forward tests: passed")


if __name__ == "__main__":
    main()
