#!/usr/bin/env python3
from __future__ import annotations

import random
import tempfile
import unittest
from pathlib import Path

import low_abundance_rescue as lr


class LowAbundanceRescueTests(unittest.TestCase):
    def test_nonbranching_paths_do_not_cross_ambiguous_junction(self):
        with tempfile.TemporaryDirectory() as td:
            gfa = Path(td) / "assembly.gfa"
            seq = "A" * 80
            gfa.write_text(
                "S\tA\t" + seq + "\tKC:f:3\n"
                "S\tB\t" + "C" * 80 + "\tKC:f:3\n"
                "S\tC\t" + "G" * 80 + "\tKC:f:3\n"
                "S\tD\t" + "T" * 80 + "\tKC:f:3\n"
                "L\tA\t+\tB\t+\t21M\tDR:i:4\n"
                "L\tB\t+\tC\t+\t21M\tDR:i:4\n"
                "L\tB\t+\tD\t+\t21M\tDR:i:4\n"
            )
            graph = lr.parse_gfa(gfa)
            paths = lr.maximal_nonbranching_paths(graph)
            self.assertFalse(any("C" in p and "D" in p for p in paths))
            self.assertFalse(any(p[:3] == ["A", "B", "C"] for p in paths))

    def test_isolated_node_is_recovered(self):
        graph = lr.GraphData(
            "x",
            21,
            {"iso": "ACGT" * 60},
            {"iso": 2.0},
            {"iso": []},
            {"iso": []},
            {},
        )
        self.assertEqual(lr.maximal_nonbranching_paths(graph), [["iso"]])

    def test_representation_fraction_separates_novel_read(self):
        k = 5
        represented = set(lr.kmers("ACGTACGTACGT", k))
        self.assertGreater(lr.represented_fraction("ACGTACGT", represented, k, 1), 0.9)
        self.assertEqual(lr.represented_fraction("TTTTTTTT", represented, k, 1), 0.0)

    def test_graph_selection_rejects_incoherent_coverage(self):
        rng = random.Random(7)
        good_seq = "".join(rng.choice("ACGT") for _ in range(320))
        bad_seq = "".join(rng.choice("ACGT") for _ in range(320))
        good = lr.GraphCandidate(
            "g", ["a", "b"], good_seq, 2.5, 2.0, 3.0, 1.5, 100, 0.8, 0, 5, 100.0
        )
        bad = lr.GraphCandidate(
            "g", ["x", "y"], bad_seq, 5.5, 1.0, 10.0, 10.0, 100, 0.8, 0, 5, 200.0
        )
        selected = lr.select_graph_candidates(
            [good, bad],
            set(),
            min_novel_kmers=20,
            min_novel_fraction=0.5,
            max_coverage_ratio=4.0,
            max_total_bases=10000,
            per_cluster_fraction=1.0,
            novel_k=31,
        )
        self.assertEqual([x.nodes for x in selected], [["a", "b"]])


if __name__ == "__main__":
    unittest.main()
