#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter

import stage23_guided_singleton_rescue as s23


def synthetic_case():
    k = 5
    seq = "ACGTTGCATGACCTGATCGTAGCTAGGCTA"
    record = s23.PriorRecord("u1", seq)
    nodes = s23.ordered_keys(seq, k)
    edges = s23.ordered_keys(seq, k + 1)
    observations = Counter({key: 2 for key in nodes})
    fragments = Counter({key: 2 for key in nodes})
    quality_sum = Counter({key: 2 * k * 35 for key in nodes})
    supporters = {}
    edge_fragments = Counter({key: 1 for key in edges})
    # Turn two disjoint internal nodes into high-quality singleton nodes with
    # different physical-fragment supporters.  Solid nodes remain on both sides.
    for pos, pair_id in ((3, 101), (8, 202)):
        key = nodes[pos]
        observations[key] = 1
        fragments[key] = 1
        quality_sum[key] = k * 35
        supporters[key] = pair_id
    return k, record, nodes, edges, observations, fragments, quality_sum, supporters, edge_fragments


def test_selects_raw_edge_constrained_multi_fragment_singletons() -> None:
    k, record, nodes, _edges, observations, fragments, quality_sum, supporters, edge_fragments = synthetic_case()
    segments, rescued, stats = s23.select_guided_segments(
        [record],
        k=k,
        observations=observations,
        fragment_count=fragments,
        quality_sum=quality_sum,
        singleton_supporter=supporters,
        edge_fragments=edge_fragments,
        singleton_quality=30.0,
        min_distinct_singleton_fragments=2,
    )
    assert len(segments) == 1
    assert nodes[3] in rescued and nodes[8] in rescued
    assert stats["selected_segments"] == 1


def test_rejects_when_a_real_edge_is_missing() -> None:
    k, record, _nodes, edges, observations, fragments, quality_sum, supporters, edge_fragments = synthetic_case()
    edge_fragments[edges[5]] = 0
    segments, _rescued, _stats = s23.select_guided_segments(
        [record],
        k=k,
        observations=observations,
        fragment_count=fragments,
        quality_sum=quality_sum,
        singleton_supporter=supporters,
        edge_fragments=edge_fragments,
        singleton_quality=30.0,
        min_distinct_singleton_fragments=2,
    )
    assert segments == []


def test_rejects_single_fragment_error_like_path() -> None:
    k, record, nodes, _edges, observations, fragments, quality_sum, supporters, edge_fragments = synthetic_case()
    supporters[nodes[8]] = 101
    segments, _rescued, stats = s23.select_guided_segments(
        [record],
        k=k,
        observations=observations,
        fragment_count=fragments,
        quality_sum=quality_sum,
        singleton_supporter=supporters,
        edge_fragments=edge_fragments,
        singleton_quality=30.0,
        min_distinct_singleton_fragments=2,
    )
    assert segments == []
    assert stats["rejected_single_fragment"] >= 1


def test_n_breaks_ordered_path() -> None:
    assert s23.ordered_keys("AAAAANCCCCC", 5) == []


def main() -> None:
    test_selects_raw_edge_constrained_multi_fragment_singletons()
    test_rejects_when_a_real_edge_is_missing()
    test_rejects_single_fragment_error_like_path()
    test_n_breaks_ordered_path()
    print("stage23 guided singleton tests: 4 passed")


if __name__ == "__main__":
    main()
