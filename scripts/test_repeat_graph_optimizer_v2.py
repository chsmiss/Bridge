#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
import unittest

import graph_path_phaser as gp
import repeat_graph_optimizer_v2 as rv2


class RepeatSeededTraversalTest(unittest.TestCase):
    def test_unique_flank_seed_carries_context_through_repeat(self):
        graph = gp.Graph()
        graph.k = 31
        graph.seqs = {uid: base * 62 for uid, base in {
            "a": "A", "b": "C", "r": "G", "c": "T", "d": "A"
        }.items()}
        graph.coverage = {"a": 10.0, "b": 8.0, "r": 40.0, "c": 10.0, "d": 8.0}
        graph.rev = {uid: uid for uid in graph.seqs}
        for uid in graph.seqs:
            graph.out[uid] = []
            graph.inc[uid] = []
        for src, dst in (("a", "r"), ("b", "r"), ("r", "c"), ("r", "d")):
            graph.out[src].append(dst)
            graph.inc[dst].append(src)
            graph.edge[(src, dst)] = gp.EdgeEvidence(5, 0, 0)
        raw = Counter({("a", "r", "c"): 10, ("b", "r", "d"): 4})
        paths, stats = rv2.resolve_context_seeded_paths(
            graph,
            raw,
            Counter(),
            Counter(),
            Counter(),
            0.70,
            4,
            50,
        )
        self.assertTrue(any(path == ["a", "r", "c"] for path in paths))
        self.assertGreater(stats["phased_extensions"], 0)


if __name__ == "__main__":
    unittest.main()
