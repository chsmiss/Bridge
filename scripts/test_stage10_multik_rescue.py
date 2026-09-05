#!/usr/bin/env python3
from __future__ import annotations

import random
import unittest

import stage10_multik_rescue as stage10


def random_seq(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class Stage10Tests(unittest.TestCase):
    def test_cross_k_consensus_is_detected(self):
        shared = random_seq(7, 420)
        raw = [
            (17, "k17", "a", shared),
            (21, "k21", "b", shared),
            (31, "k31", "c", random_seq(8, 420)),
        ]
        cs = stage10.annotate_multik_candidates(raw, set(), set())
        pair = [c for c in cs if c.k in (17, 21)]
        self.assertTrue(all(c.cross_k_sources >= 1 for c in pair))
        self.assertTrue(all(c.cross_k_fraction > 0.9 for c in pair))
        lone = next(c for c in cs if c.k == 31)
        self.assertEqual(lone.cross_k_sources, 0)

    def test_strict_rejects_single_k_even_if_novel(self):
        seq = random_seq(9, 430)
        cs = stage10.annotate_multik_candidates([(31, "k31", "x", seq)], set(), set())
        selected = stage10.select_multik_candidates(
            cs,
            set(),
            min_novel_kmers=64,
            min_novel_fraction=0.7,
            min_cross_sources=1,
            min_cross_fraction=0.3,
            max_total_bases=10000,
            max_fraction_per_k=1.0,
            allow_strong_single_k=False,
        )
        self.assertEqual(selected, [])

    def test_balanced_can_keep_very_strong_single_k(self):
        seq = random_seq(10, 430)
        cs = stage10.annotate_multik_candidates([(25, "k25", "x", seq)], set(), set())
        selected = stage10.select_multik_candidates(
            cs,
            set(),
            min_novel_kmers=40,
            min_novel_fraction=0.5,
            min_cross_sources=1,
            min_cross_fraction=0.2,
            max_total_bases=10000,
            max_fraction_per_k=1.0,
            allow_strong_single_k=True,
        )
        self.assertEqual(len(selected), 1)

    def test_fresh_novel_filter_deduplicates_cross_k_copies(self):
        seq = random_seq(11, 450)
        raw = [(17, "k17", "a", seq), (21, "k21", "b", seq)]
        cs = stage10.annotate_multik_candidates(raw, set(), set())
        selected = stage10.select_multik_candidates(
            cs,
            set(),
            min_novel_kmers=40,
            min_novel_fraction=0.5,
            min_cross_sources=1,
            min_cross_fraction=0.2,
            max_total_bases=10000,
            max_fraction_per_k=1.0,
            allow_strong_single_k=False,
        )
        self.assertEqual(len(selected), 1)


if __name__ == "__main__":
    unittest.main()
