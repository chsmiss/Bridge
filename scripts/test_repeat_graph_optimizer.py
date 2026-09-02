#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
import unittest

import graph_path_phaser as gp
import repeat_graph_optimizer as rg


class RepeatGraphOptimizerTest(unittest.TestCase):
    def graph(self) -> gp.Graph:
        graph = gp.Graph()
        graph.k = 31
        graph.seqs = {
            "u0": "A" * 62,
            "u1": "C" * 62,
            "u2": "G" * 62,
            "u3": "T" * 62,
        }
        graph.coverage = {"u0": 10.0, "u1": 10.0, "u2": 1.0, "u3": 10.0}
        graph.rev = {uid: uid for uid in graph.seqs}
        graph.out["u0"] = ["u1", "u2"]
        graph.out["u1"] = ["u3"]
        graph.inc["u1"] = ["u0"]
        graph.inc["u2"] = ["u0"]
        graph.inc["u3"] = ["u1"]
        graph.edge[("u0", "u1")] = gp.EdgeEvidence(10, 0, 1)
        graph.edge[("u0", "u2")] = gp.EdgeEvidence(1, 0, 0)
        graph.edge[("u1", "u3")] = gp.EdgeEvidence(10, 0, 1)
        return graph

    def test_simplification_masks_unsupported_short_tip(self):
        graph = self.graph()
        context = Counter({("u0", "u1", "u3"): 4})
        simplified, stats = rg.simplify_graph(graph, context)
        self.assertIn("u1", simplified.out["u0"])
        self.assertNotIn("u2", simplified.out["u0"])
        self.assertGreaterEqual(stats["masked_edges"], 1)

    def test_long_context_vetoes_pruning(self):
        graph = self.graph()
        context = Counter({("u0", "u1", "u3"): 4, ("u3", "u0", "u2"): 2})
        simplified, _ = rg.simplify_graph(graph, context)
        self.assertIn("u2", simplified.out["u0"])

    def test_join_paths_uses_only_existing_unique_bridge(self):
        graph = self.graph()
        merged = rg.join_paths(["u0"], ["u3"], graph, max_edges=3, max_span=200)
        self.assertEqual(merged, ["u0", "u1", "u3"])


if __name__ == "__main__":
    unittest.main()
