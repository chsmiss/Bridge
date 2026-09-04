#!/usr/bin/env python3
from __future__ import annotations

import tempfile
from pathlib import Path

import graph_path_phaser as gp
import stage16_root_cause as s16
import stage17_frontier_recovery as s17


def write_fastq(path: Path, seqs: list[str]) -> None:
    with path.open("w") as handle:
        for i, seq in enumerate(seqs):
            handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def make_graph() -> gp.Graph:
    g = gp.Graph()
    g.k = 5
    g.seqs = {
        "A": "ACGTACCCG",
        "B": "CCCGTTTAA",
        "X": "GGGGGAAAA",
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


def test_hybrid_keeps_all_exact_contexts() -> None:
    graph = make_graph()
    exact = gp.KmerIndex(graph, 5)
    fallback = s16.SparseGraphIndex(graph, 3, max_occurrences=4)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        r1 = root / "r1.fastq"
        r2 = root / "r2.fastq"
        seq = "ACGTACCCGTTTAA"
        write_fastq(r1, [seq])
        write_fastq(r2, [seq])
        exact_ctx, fallback_ctx, hybrid, stats = s17.collect_hybrid_contexts(
            graph,
            exact,
            [("k3", fallback, 2, 4)],
            r1,
            r2,
            max_context=6,
        )
        assert exact_ctx
        assert not fallback_ctx
        for key, value in exact_ctx.items():
            assert hybrid[key] >= value
        assert stats["exact_threaded_reads"] == 2
        assert stats["k3_threaded_reads"] == 0


def test_frontier_recruits_mate_borne_extension() -> None:
    source = "ACGTTGCAACGTAGCTTGCAACGTTGCAACGTAGCTTGCA"
    pairs = [(source, source[::-1])]
    assignments = {0: 0}
    frontier = s17.locus_unique_frontier(
        pairs,
        assignments,
        set(),
        k=19,
        stride=3,
        min_count=1,
    )
    assert len(frontier) >= 2
    mers = list(frontier)[:2]
    recruited = mers[0] + "AA" + mers[1]
    pairs.append((recruited, "T" * len(recruited)))
    added, ambiguous = s17.recruit_frontier_round(
        pairs,
        assignments,
        frontier,
        k=19,
        stride=3,
        min_hits=2,
        margin=1,
    )
    assert ambiguous == 0
    assert added == 1
    assert assignments[1] == 0


def test_frontier_drops_cross_locus_kmers() -> None:
    shared = "ACGTTGCAACGTAGCTTGCAACGTTGCAACGTAGCTTGCA"
    pairs = [(shared, shared), (shared, shared)]
    assignments = {0: 0, 1: 1}
    frontier = s17.locus_unique_frontier(
        pairs,
        assignments,
        set(),
        k=19,
        stride=3,
        min_count=1,
    )
    assert not frontier


def test_path_support_prefers_longer_context() -> None:
    from collections import Counter

    path = ["A", "B", "C", "D"]
    ctx = Counter({("A", "B"): 5, ("A", "B", "C"): 1})
    support, span = s17.path_context_support(path, ctx)
    assert span == 3
    assert support == 1


def main() -> None:
    tests = [
        test_hybrid_keeps_all_exact_contexts,
        test_frontier_recruits_mate_borne_extension,
        test_frontier_drops_cross_locus_kmers,
        test_path_support_prefers_longer_context,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)


if __name__ == "__main__":
    main()
