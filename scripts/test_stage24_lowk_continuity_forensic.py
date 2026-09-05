#!/usr/bin/env python3
from __future__ import annotations

import tempfile
from pathlib import Path

import stage24_lowk_continuity_forensic as s24


def test_consecutive_true_runs() -> None:
    assert s24.consecutive_true_runs([False, True, True, False, True]) == [(1, 2), (4, 4)]
    assert s24.consecutive_true_runs([True, True]) == [(0, 1)]
    assert s24.consecutive_true_runs([False, False]) == []


def test_exact_nonbranching_spelling() -> None:
    seq = "ACGTACGT"
    k = 3
    keys = list(s24.oriented_keys(seq, k))
    nodes = set(keys)
    edges = set(zip(keys, keys[1:]))
    out = {}
    inc = {}
    for left, right in edges:
        out.setdefault(left, set()).add(right)
        inc.setdefault(right, set()).add(left)
    metrics = s24.evaluate_k21_spelling(seq, nodes, edges, out, inc, k=k)
    assert metrics["k21_exact_path"] is True
    assert metrics["k21_nonbranching"] is True
    assert s24.classify_run(metrics, True) == "exact_nonbranching_both_k31_anchors"


def test_branching_spelling_is_not_direct_actionable() -> None:
    seq = "ACGTACGT"
    k = 3
    keys = list(s24.oriented_keys(seq, k))
    nodes = set(keys)
    edges = set(zip(keys, keys[1:]))
    out = {}
    inc = {}
    for left, right in edges:
        out.setdefault(left, set()).add(right)
        inc.setdefault(right, set()).add(left)
    internal = keys[len(keys) // 2]
    out.setdefault(internal, set()).add(999999)
    metrics = s24.evaluate_k21_spelling(seq, nodes, edges, out, inc, k=k)
    assert metrics["k21_exact_path"] is True
    assert metrics["k21_nonbranching"] is False
    assert s24.classify_run(metrics, True) == "exact_branching_both_k31_anchors"


def test_missing_edge_separates_graph_break_from_sequence_loss() -> None:
    seq = "ACGTACGT"
    k = 3
    keys = list(s24.oriented_keys(seq, k))
    all_edges = list(zip(keys, keys[1:]))
    nodes = set(keys)
    edges = set(all_edges[:-1])
    out = {}
    inc = {}
    for left, right in edges:
        out.setdefault(left, set()).add(right)
        inc.setdefault(right, set()).add(left)
    metrics = s24.evaluate_k21_spelling(seq, nodes, edges, out, inc, k=k)
    assert metrics["k21_exact_path"] is False
    assert s24.classify_run(metrics, True) == "k21_graph_broken"

    sparse_nodes = {keys[0]}
    sparse = s24.evaluate_k21_spelling(seq, sparse_nodes, set(), {}, {}, k=k)
    assert s24.classify_run(sparse, True) == "k21_sequence_missing"


def test_graph_topology_reads_directed_unitig_transitions() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        gfa = Path(tmp) / "graph.gfa"
        gfa.write_text("H\tVN:Z:1.0\nS\tu1\tACGTAC\nS\tu2\tGTACGA\n")
        nodes, edges, out, inc, records = s24.graph_topology(gfa, k=3)
        expected = list(s24.oriented_keys("ACGTAC", 3))
        assert records == 2
        assert set(expected).issubset(nodes)
        assert all(pair in edges for pair in zip(expected, expected[1:]))
        assert expected[1] in out[expected[0]]
        assert expected[0] in inc[expected[1]]


def main() -> None:
    test_consecutive_true_runs()
    test_exact_nonbranching_spelling()
    test_branching_spelling_is_not_direct_actionable()
    test_missing_edge_separates_graph_break_from_sequence_loss()
    test_graph_topology_reads_directed_unitig_transitions()
    print("stage24 low-k forensic tests passed")


if __name__ == "__main__":
    main()
