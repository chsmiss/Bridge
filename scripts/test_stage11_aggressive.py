#!/usr/bin/env python3
from __future__ import annotations

import unittest
from collections import Counter

import fill_scaffold_gaps_multik as fm
import graph_path_phaser as gp
import stage11_aggressive_rescue as s11


class Stage11AggressiveTests(unittest.TestCase):
    def test_multik_gap_fill_requires_consensus(self):
        fill, ks = fm.consensus_fill({17: "ACGT", 21: "ACGT", 25: "ACGA"}, 2)
        self.assertEqual(fill, "ACGT")
        self.assertEqual(ks, [17, 21])
        fill2, ks2 = fm.consensus_fill({17: "ACGT", 21: "ACGA", 25: None}, 2)
        self.assertIsNone(fill2)
        self.assertEqual(ks2, [])

    def test_unresolved_n_is_split_not_emitted_as_scaffold(self):
        records = [("x", "A" * 220 + "N" * 50 + "C" * 240)]
        out = fm.split_unresolved(records, 200)
        self.assertEqual([len(seq) for _name, seq in out], [220, 240])
        self.assertTrue(all("N" not in seq for _name, seq in out))

    def make_graph(self) -> gp.Graph:
        graph = gp.Graph()
        graph.seqs = {
            "X": "A" * 80,
            "A": "C" * 80,
            "B": "G" * 80,
            "C": "T" * 80,
        }
        graph.coverage = {"X": 10.0, "A": 10.0, "B": 9.0, "C": 2.0}
        graph.out = {"X": ["A"], "A": ["B", "C"], "B": [], "C": []}
        graph.inc = {"X": [], "A": ["X"], "B": ["A"], "C": ["A"]}
        graph.rev = {uid: uid for uid in graph.seqs}
        graph.edge = {
            ("X", "A"): gp.EdgeEvidence(4, 0, 0),
            ("A", "B"): gp.EdgeEvidence(2, 0, 2),
            ("A", "C"): gp.EdgeEvidence(2, 0, 0),
        }
        return graph

    def test_abundance_flow_ranks_supported_dominant_branch(self):
        graph = self.make_graph()
        raw = Counter({("X", "A", "B"): 2})
        ranked = s11.rank_flow_candidates(
            graph,
            ["X", "A"],
            ["B", "C"],
            set(),
            True,
            raw,
            Counter(),
            Counter(),
            Counter(),
        )
        self.assertTrue(ranked)
        self.assertEqual(ranked[0].choice.uid, "B")
        self.assertGreater(ranked[0].sibling_share, 0.8)
        self.assertGreaterEqual(ranked[0].evidence_channels, 2)

    def test_coverage_alone_cannot_force_branch(self):
        graph = self.make_graph()
        graph.edge[("A", "B")] = gp.EdgeEvidence(0, 0, 0)
        graph.edge[("A", "C")] = gp.EdgeEvidence(0, 0, 0)
        ranked = s11.rank_flow_candidates(
            graph,
            ["X", "A"],
            ["B", "C"],
            set(),
            True,
            Counter(),
            Counter(),
            Counter(),
            Counter(),
        )
        self.assertEqual(ranked, [])


if __name__ == "__main__":
    unittest.main()
