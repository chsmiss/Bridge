#!/usr/bin/env python3
from __future__ import annotations

import stage17_breakpoint_oracle as oracle


def test_canonical_kmers_are_strand_invariant() -> None:
    seq = "ACGTTGCAACGTAGCTTGCAACGTTGCAACGTAGCTTGCA"
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert set(oracle.rolling_canonical_keys(seq, 21)) == set(
        oracle.rolling_canonical_keys(rc, 21)
    )


def test_uncovered_intervals() -> None:
    aligned = oracle.merge_intervals([(100, 300), (320, 500), (800, 900)], join_gap=50)
    assert aligned == [(100, 500), (800, 900)]
    assert oracle.uncovered_intervals(1200, aligned, min_gap=200) == [(500, 800), (900, 1200)]


def test_classification_order() -> None:
    assert oracle.classify(0.05, 0.02, 0.0, 0.0, 0.0) == "sampling_or_no_read_support"
    assert oracle.classify(0.8, 0.7, 0.1, 0.0, 0.0) == "low_k_graph_loss"
    assert oracle.classify(0.8, 0.7, 0.8, 0.2, 0.0) == "k21_to_k31_propagation_loss"
    assert oracle.classify(0.8, 0.7, 0.8, 0.8, 0.2) == "emission_or_path_selection_loss"


def main() -> None:
    tests = [
        test_canonical_kmers_are_strand_invariant,
        test_uncovered_intervals,
        test_classification_order,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)


if __name__ == "__main__":
    main()
