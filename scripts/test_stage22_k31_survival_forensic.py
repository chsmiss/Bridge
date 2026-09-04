#!/usr/bin/env python3
from __future__ import annotations

import stage22_k31_survival_forensic as forensic


def test_canonical_reverse_complement() -> None:
    left = forensic.canonical_key("A" * 31 + "C")
    right = forensic.canonical_key("G" + "T" * 31)
    assert left == right


def test_sampled_intervals() -> None:
    assert forensic.sampled_intervals(0, 1000) == [(0, 1000)]
    sampled = forensic.sampled_intervals(0, 5000)
    assert sampled[0] == (0, 500)
    assert sampled[-1] == (4500, 5000)
    assert len(sampled) == 3


def test_rolling_keys_break_on_n() -> None:
    observed = list(forensic.rolling_keys("A" * 31 + "N" + "C" * 31, 31))
    assert len(observed) == 2
    assert observed[0][0] == 0
    assert observed[1][0] == 32


def main() -> None:
    test_canonical_reverse_complement()
    test_sampled_intervals()
    test_rolling_keys_break_on_n()
    print("stage22 forensic tests: 3 passed")


if __name__ == "__main__":
    main()
