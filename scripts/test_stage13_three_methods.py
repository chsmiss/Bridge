#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter

import stage13_three_methods as s13


def make_graph():
    graph = s13.gp.Graph()
    for uid, cov in [("S", 10.0), ("A", 8.0), ("B", 2.0), ("T", 10.0)]:
        graph.seqs[uid] = "A" * 80
        graph.coverage[uid] = cov
        graph.rev[uid] = uid
    graph.out["S"] = ["A", "B"]
    graph.out["A"] = ["T"]
    graph.out["B"] = ["T"]
    graph.out["T"] = []
    graph.inc["S"] = []
    graph.inc["A"] = ["S"]
    graph.inc["B"] = ["S"]
    graph.inc["T"] = ["A", "B"]
    return graph


def test_nnls_preserves_low_abundance_path():
    paths = [["S", "A", "T"], ["S", "B", "T"]]
    values, rmse = s13.projected_nnls(
        paths, {"A": 8.0, "B": 2.0}, {"A": 100, "B": 100}, l1=0.0
    )
    assert abs(values[0] - 8.0) < 1e-6
    assert abs(values[1] - 2.0) < 1e-6
    assert rmse < 1e-6


def test_bubble_discovery_finds_two_paths():
    graph = make_graph()
    bubbles = s13.discover_bubbles(
        graph, max_depth=4, max_paths=4, max_branch=4, max_bubbles=10
    )
    assert len(bubbles) == 1
    bubble = bubbles[0]
    assert bubble.source == "S"
    assert bubble.sink == "T"
    assert {tuple(path) for path in bubble.paths} == {
        ("S", "A", "T"),
        ("S", "B", "T"),
    }


def test_thread_support_accepts_partial_path_context():
    counter = Counter({("S", "A", "X"): 3, ("A", "X"): 5})
    assert s13.path_thread_support(counter, ["S", "A", "X", "T"]) == 5


def test_physical_evidence_counts_supported_edges():
    graph = make_graph()
    graph.edge[("S", "A")] = s13.gp.EdgeEvidence(direct=2, gapped=0, pairs=1)
    graph.edge[("A", "T")] = s13.gp.EdgeEvidence()
    physical, edges, pairs = s13.path_physical_evidence(graph, ["S", "A", "T"])
    assert (physical, edges, pairs) == (1, 2, 1)


def main():
    tests = [
        test_nnls_preserves_low_abundance_path,
        test_bubble_discovery_finds_two_paths,
        test_thread_support_accepts_partial_path_context,
        test_physical_evidence_counts_supported_edges,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
