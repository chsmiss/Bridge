#!/usr/bin/env python3
"""Weighted-routing variant of BridgeAsm component-local k-lifting.

V3 required a routing seed k-mer to be unique to exactly one graph neighborhood.
That is too strict precisely around repeats/strain bubbles, where useful anchors are
often shared by a small number of nearby neighborhoods. V4 keeps bounded shared
signatures and divides each seed's vote across its memberships. Ambiguous read
pairs may then be copied to the two strongest neighborhoods, preserving the soft
assignment semantics without letting ubiquitous k-mers dominate routing.
"""
from __future__ import annotations

import argparse
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO

import adaptive_k_local_v2 as v2
import local_k_lift_v3 as v3


def build_weighted_signature_index(
    neighborhoods: list[v2.Neighborhood],
    segments: dict[str, v2.Segment],
    seed_k: int,
    stride: int,
    max_memberships: int = 4,
) -> dict[str, tuple[int, ...]]:
    memberships: dict[str, set[int]] = defaultdict(set)
    for neighborhood in neighborhoods:
        local: set[str] = set()
        for node in neighborhood.nodes:
            local.update(v2.canonical_kmers(segments[node].sequence, seed_k, stride))
        for kmer in local:
            memberships[kmer].add(neighborhood.identifier)
    return {
        kmer: tuple(sorted(ids))
        for kmer, ids in memberships.items()
        if 1 <= len(ids) <= max_memberships
    }


def route_pairs_weighted_soft(
    read1: Path,
    read2: Path,
    output_dir: Path,
    neighborhoods: list[v2.Neighborhood],
    signatures: dict[str, tuple[int, ...]],
    seed_k: int,
    stride: int,
    min_hits: int,
    min_top_fraction: float,
    soft_second_fraction: float,
) -> tuple[int, int, dict[int, int]]:
    handles: dict[int, tuple[TextIO, TextIO]] = {}
    soft_counts: Counter[int] = Counter()
    by_id = {item.identifier: item for item in neighborhoods}
    for item in neighborhoods:
        directory = output_dir / f"neighborhood_{item.identifier:03d}" / "reads"
        directory.mkdir(parents=True, exist_ok=True)
        handles[item.identifier] = (
            gzip.open(directory / "R1.fastq.gz", "wt"),
            gzip.open(directory / "R2.fastq.gz", "wt"),
        )

    routed_pairs = 0
    assignments = 0
    left_iter = v2.fastq_records(read1)
    right_iter = v2.fastq_records(read2)
    while True:
        try:
            left = next(left_iter)
        except StopIteration:
            left = None
        try:
            right = next(right_iter)
        except StopIteration:
            right = None
        if left is None and right is None:
            break
        if left is None or right is None:
            raise ValueError("paired FASTQ files contain different record counts")

        scores: Counter[int] = Counter()
        total_evidence = 0.0
        for sequence in (left[1], right[1]):
            for kmer in v2.canonical_kmers(sequence, seed_k, stride):
                identifiers = signatures.get(kmer)
                if not identifiers:
                    continue
                weight = 1.0 / len(identifiers)
                for identifier in identifiers:
                    scores[identifier] += weight
                total_evidence += 1.0

        ranked = scores.most_common(2)
        if not ranked or ranked[0][1] < float(min_hits) or total_evidence == 0.0:
            continue
        if ranked[0][1] / total_evidence < min_top_fraction:
            continue

        selected = [ranked[0][0]]
        if (
            len(ranked) > 1
            and ranked[1][1] >= float(min_hits)
            and ranked[1][1] >= ranked[0][1] * soft_second_fraction
        ):
            selected.append(ranked[1][0])

        routed_pairs += 1
        for identifier in selected:
            out1, out2 = handles[identifier]
            v2.write_record(out1, left)
            v2.write_record(out2, right)
            by_id[identifier].pair_count += 1
            if len(selected) > 1:
                soft_counts[identifier] += 1
            assignments += 1

    for out1, out2 in handles.values():
        out1.close()
        out2.close()
    return routed_pairs, assignments, dict(soft_counts)


def main() -> None:
    # V3's main owns the assembly/promotion/reconciliation pipeline. Patch only
    # the two routing hooks so this experiment changes one variable at a time.
    original_builder = v3.v2.build_signature_index
    original_router = v3.route_pairs_soft

    def builder(
        neighborhoods: list[v2.Neighborhood],
        segments: dict[str, v2.Segment],
        seed_k: int,
        stride: int,
    ) -> dict[str, tuple[int, ...]]:
        return build_weighted_signature_index(
            neighborhoods, segments, seed_k, stride, max_memberships=4
        )

    try:
        v3.v2.build_signature_index = builder
        v3.route_pairs_soft = route_pairs_weighted_soft
        v3.main()
    finally:
        v3.v2.build_signature_index = original_builder
        v3.route_pairs_soft = original_router


if __name__ == "__main__":
    main()
