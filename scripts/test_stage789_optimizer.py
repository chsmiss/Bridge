#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
import unittest

import graph_path_phaser as gp
import stage789_optimizer as s789


class Stage789LookaheadTest(unittest.TestCase):
    @staticmethod
    def make_graph() -> gp.Graph:
        graph = gp.Graph()
        graph.k = 31
        graph.seqs = {uid: ("ACGT" * 40) for uid in ("a", "r", "x", "y", "c", "d")}
        graph.coverage = {uid: 10.0 for uid in graph.seqs}
        graph.rev = {uid: uid for uid in graph.seqs}
        for uid in graph.seqs:
            graph.out[uid] = []
            graph.inc[uid] = []
        for src, dst, direct in (
            ("a", "r", 5),
            ("r", "x", 2),
            ("r", "y", 2),
            ("x", "c", 1),
            ("y", "d", 1),
        ):
            graph.out[src].append(dst)
            graph.inc[dst].append(src)
            graph.edge[(src, dst)] = gp.EdgeEvidence(direct, 0, 0)
        return graph

    def choose(self, raw: Counter[tuple[str, ...]]):
        return s789.choose_extension_lookahead(
            self.make_graph(),
            ["a", "r"],
            ["x", "y"],
            {"a", "r"},
            True,
            raw,
            Counter(),
            Counter(),
            Counter(),
            0.70,
            4,
            3,
            4,
            0.70,
            0.60,
            1.15,
        )

    def test_downstream_context_rescues_stalled_branch(self):
        choice, rescued = self.choose(Counter({("r", "x", "c"): 8, ("r", "y", "d"): 1}))
        self.assertIsNotNone(choice)
        assert choice is not None
        self.assertEqual(choice.uid, "x")
        self.assertTrue(rescued)

    def test_symmetric_downstream_context_stays_unresolved(self):
        choice, rescued = self.choose(Counter({("r", "x", "c"): 5, ("r", "y", "d"): 5}))
        self.assertIsNone(choice)
        self.assertFalse(rescued)

    def test_weak_root_edge_is_not_rescued(self):
        graph = self.make_graph()
        graph.edge[("r", "x")] = gp.EdgeEvidence(0, 0, 0)
        choice, rescued = s789.choose_extension_lookahead(
            graph,
            ["a", "r"],
            ["x", "y"],
            {"a", "r"},
            True,
            Counter({("r", "x", "c"): 20, ("r", "y", "d"): 1}),
            Counter(), Counter(), Counter(),
            0.70, 4, 3, 4, 0.70, 0.60, 1.15,
        )
        self.assertIsNone(choice)
        self.assertFalse(rescued)


if __name__ == "__main__":
    unittest.main()
