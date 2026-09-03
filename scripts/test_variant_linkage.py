#!/usr/bin/env python3
import unittest
from collections import Counter

import graph_path_phaser as gp
import variant_linkage as vl


def toy_graph():
    g = gp.Graph()
    for uid in ["A", "B", "C", "P", "Q"]:
        g.seqs[uid] = "A" * 80
        g.rev[uid] = uid
    g.out["A"] = ["B", "C"]
    g.inc["B"] = ["A"]
    g.inc["C"] = ["A"]
    return g


class VariantLinkageTests(unittest.TestCase):
    def test_history_excludes_current_branch_source(self):
        self.assertEqual(vl.history_context(["P", "A"], True, 3), ["P"])
        self.assertEqual(vl.history_context(["A", "P", "Q"], False, 2), ["P", "Q"])

    def test_dominant_haplotype_scores_correct_child(self):
        g = toy_graph()
        linkage = vl.MarkerLinkage(
            Counter({("B", "P"): 8, ("C", "P"): 1}), Counter(), 10, 9, 21, 2
        )
        scores = vl.score_candidates(g, linkage, ["P", "A"], ["B", "C"], True)
        self.assertEqual(scores[0].uid, "B")
        self.assertGreater(scores[0].share, 0.85)

    def test_guard_vetoes_contradicted_flow_choice(self):
        g = toy_graph()
        linkage = vl.MarkerLinkage(
            Counter({("B", "P"): 9, ("C", "P"): 1}), Counter(), 10, 9, 21, 2
        )
        mode, uid, _ = vl.linkage_decision(
            g, linkage, ["P", "A"], ["B", "C"], True, "C"
        )
        self.assertEqual(mode, "veto")
        self.assertIsNone(uid)

    def test_guard_rescues_strong_unresolved_branch(self):
        g = toy_graph()
        linkage = vl.MarkerLinkage(
            Counter({("B", "P"): 7, ("C", "P"): 2}), Counter(), 10, 9, 21, 2
        )
        mode, uid, _ = vl.linkage_decision(
            g, linkage, ["P", "A"], ["B", "C"], True, None
        )
        self.assertEqual(mode, "rescue")
        self.assertEqual(uid, "B")

    def test_weak_linkage_does_not_override(self):
        g = toy_graph()
        linkage = vl.MarkerLinkage(Counter({("B", "P"): 1}), Counter(), 10, 2, 21, 2)
        mode, uid, _ = vl.linkage_decision(
            g, linkage, ["P", "A"], ["B", "C"], True, None
        )
        self.assertEqual(mode, "none")
        self.assertIsNone(uid)


if __name__ == "__main__":
    unittest.main()
