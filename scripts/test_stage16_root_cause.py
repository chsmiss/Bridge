#!/usr/bin/env python3
from __future__ import annotations

import tempfile
from pathlib import Path

import graph_path_phaser as gp
import stage16_root_cause as s16


def make_graph() -> gp.Graph:
    g = gp.Graph()
    g.k = 5
    g.seqs = {
        "A": "ACGTACCCG",
        "B": "CCCGTTTAA",
        "X": "ACGTAGGGC",
        "Z": "TTTTCCCCA",
    }
    g.coverage = {uid: 1.0 for uid in g.seqs}
    for uid in g.seqs:
        g.out[uid] = []
        g.inc[uid] = []
        g.rev[uid] = uid
    g.out["A"] = ["B"]
    g.inc["B"] = ["A"]
    g.edge[("A", "B")] = gp.EdgeEvidence(2, 0, 1)
    return g


def test_ambiguous_anchors_resolve_by_topology() -> None:
    graph = make_graph()
    events = [
        s16.AnchorEvent(4, ("A", "X")),
        s16.AnchorEvent(35, ("B", "Z")),
    ]
    chain = s16.chain_anchor_events(events, graph)
    assert chain is not None
    assert chain.path == ("A", "B")
    assert chain.ambiguous_anchors == 2


def test_single_ambiguous_anchor_is_not_threaded() -> None:
    graph = make_graph()
    chain = s16.chain_anchor_events([s16.AnchorEvent(10, ("A", "X"))], graph)
    assert chain is None


def test_sparse_index_drops_overflow_repeats() -> None:
    graph = gp.Graph()
    graph.k = 5
    graph.seqs = {f"u{i}": "AAAAACCCCC" for i in range(4)}
    graph.coverage = {uid: 1.0 for uid in graph.seqs}
    for uid in graph.seqs:
        graph.out[uid] = []
        graph.inc[uid] = []
        graph.rev[uid] = uid
    index = s16.SparseGraphIndex(graph, 5, max_occurrences=2)
    key = next(gp.rolling_keys("AAAAA", 5))[1]
    assert key not in index.unique
    assert key not in index.ambig
    assert index.dropped_repetitive_keys > 0


def base4_word(value: int, length: int = 21) -> str:
    alphabet = "ACGT"
    chars = []
    for _ in range(length):
        chars.append(alphabet[value & 3])
        value >>= 2
    return "".join(chars)


def test_seed_signatures_are_not_truncated_to_36() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        seeds = root / "seeds.fasta"
        baseline = root / "baseline.fasta"
        baseline.write_text(">base\n" + "G" * 100 + "\n")
        with seeds.open("w") as handle:
            for i in range(48):
                mer = base4_word(1000 + i)
                handle.write(f">s{i}\n{mer}TTGCA{mer}\n")
        chosen, signatures, _sets, stats = s16.build_all_seed_signatures(
            seeds, baseline, k=21, min_signature_kmers=1
        )
        assert len(chosen) > 36
        assert stats["input_seed_records"] == 48
        assert signatures


def main() -> None:
    tests = [
        test_ambiguous_anchors_resolve_by_topology,
        test_single_ambiguous_anchor_is_not_threaded,
        test_sparse_index_drops_overflow_repeats,
        test_seed_signatures_are_not_truncated_to_36,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)


if __name__ == "__main__":
    main()
