#!/usr/bin/env python3
"""Long-range read-pair haplotype linkage for graph branch validation.

The module derives orientation-agnostic unitig-family markers from unique exact
k-mers. A physical read pair votes for a unitig family only when it contains
multiple distinct family-unique markers. Pair co-occurrence between an upstream
path family and a branch-child family becomes a long-range haplotype signal.

This is deliberately complementary to GFA DR/GR/PE tags: the current branch
source is excluded from scoring, so evidence must link the candidate to earlier
path history rather than merely restating local edge support.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import graph_path_phaser as gp


@dataclass
class MarkerLinkage:
    pair_counts: Counter[tuple[str, str]]
    family_depth: Counter[str]
    fragments: int
    informative_fragments: int
    k: int
    min_markers: int


@dataclass(frozen=True)
class LinkageScore:
    uid: str
    support: float
    share: float
    anchors: int


def family(graph: gp.Graph, uid: str) -> str:
    rev = graph.rev.get(uid, uid)
    return uid if uid <= rev else rev


def _unique_family_markers(
    sequence: str,
    graph: gp.Graph,
    index: gp.KmerIndex,
    stride: int,
) -> dict[str, set[int]]:
    markers: dict[str, set[int]] = defaultdict(set)
    for oriented in (sequence.upper(), gp.rc(sequence)):
        for pos, key in gp.rolling_keys(oriented, index.k):
            if stride > 1 and pos % stride:
                continue
            uid = index.unique.get(key)
            if uid is None:
                continue
            markers[family(graph, uid)].add(key)
    return markers


def fragment_families(
    left: str,
    right: str | None,
    graph: gp.Graph,
    index: gp.KmerIndex,
    *,
    stride: int = 2,
    min_markers: int = 2,
    max_families: int = 8,
) -> set[str]:
    markers: dict[str, set[int]] = defaultdict(set)
    for seq in (left, right):
        if not seq:
            continue
        for fam, keys in _unique_family_markers(seq, graph, index, stride).items():
            markers[fam].update(keys)
    hits = {fam for fam, keys in markers.items() if len(keys) >= min_markers}
    return hits if len(hits) <= max_families else set()


def collect_linkage(
    graph: gp.Graph,
    read1: Path,
    read2: Path,
    *,
    k: int = 21,
    stride: int = 2,
    min_markers: int = 2,
    max_families: int = 8,
) -> MarkerLinkage:
    index = gp.KmerIndex(graph, k)
    pair_counts: Counter[tuple[str, str]] = Counter()
    family_depth: Counter[str] = Counter()
    fragments = informative = 0
    for (_name1, seq1), (_name2, seq2) in zip(gp.read_fastq(read1), gp.read_fastq(read2)):
        fragments += 1
        hits = sorted(
            fragment_families(
                seq1,
                seq2,
                graph,
                index,
                stride=stride,
                min_markers=min_markers,
                max_families=max_families,
            )
        )
        if not hits:
            continue
        informative += 1
        family_depth.update(hits)
        for i, left in enumerate(hits):
            for right in hits[i + 1 :]:
                pair_counts[(left, right)] += 1
    return MarkerLinkage(
        pair_counts=pair_counts,
        family_depth=family_depth,
        fragments=fragments,
        informative_fragments=informative,
        k=k,
        min_markers=min_markers,
    )


def pair_support(linkage: MarkerLinkage, left: str, right: str) -> int:
    if left == right:
        return 0
    key = (left, right) if left < right else (right, left)
    return linkage.pair_counts.get(key, 0)


def history_context(history: list[str], forward: bool, depth: int) -> list[str]:
    if forward:
        nodes = list(reversed(history[:-1]))
    else:
        nodes = history[1:]
    out: list[str] = []
    seen: set[str] = set()
    for uid in nodes:
        if uid in seen:
            continue
        seen.add(uid)
        out.append(uid)
        if len(out) >= depth:
            break
    return out


def score_candidates(
    graph: gp.Graph,
    linkage: MarkerLinkage,
    history: list[str],
    candidates: Iterable[str],
    forward: bool,
    *,
    history_depth: int = 3,
) -> list[LinkageScore]:
    context = history_context(history, forward, history_depth)
    raw: list[tuple[str, float, int]] = []
    for uid in candidates:
        cand_family = family(graph, uid)
        support = 0.0
        anchors = 0
        seen_history_families: set[str] = set()
        for distance, hist_uid in enumerate(context, 1):
            hist_family = family(graph, hist_uid)
            if hist_family == cand_family or hist_family in seen_history_families:
                continue
            seen_history_families.add(hist_family)
            count = pair_support(linkage, hist_family, cand_family)
            if count <= 0:
                continue
            weight = 1.0 + 0.20 * min(distance - 1, 2)
            support += count * weight
            anchors += 1
        raw.append((uid, support, anchors))
    total = sum(item[1] for item in raw)
    scores = [
        LinkageScore(uid, support, support / total if total > 0 else 0.0, anchors)
        for uid, support, anchors in raw
    ]
    scores.sort(key=lambda item: (-item.support, -item.anchors, item.uid))
    return scores


def linkage_decision(
    graph: gp.Graph,
    linkage: MarkerLinkage,
    history: list[str],
    candidates: list[str],
    forward: bool,
    proposed: str | None,
    *,
    min_total: float = 3.0,
    min_rescue_support: float = 3.0,
    rescue_share: float = 0.67,
    rescue_margin: float = 1.50,
    veto_share: float = 0.25,
    veto_margin: float = 2.00,
) -> tuple[str, str | None, list[LinkageScore]]:
    scores = score_candidates(graph, linkage, history, candidates, forward)
    if not scores:
        return ("accept" if proposed is not None else "none", proposed, scores)
    by_uid = {item.uid: item for item in scores}
    total = sum(item.support for item in scores)
    best = scores[0]
    second = scores[1].support if len(scores) > 1 else 0.0

    if proposed is not None:
        current = by_uid.get(proposed, LinkageScore(proposed, 0.0, 0.0, 0))
        if total >= min_total and best.uid != proposed:
            contradicted = (
                current.share < veto_share
                and best.support >= max(min_rescue_support, current.support * veto_margin)
            )
            if contradicted:
                return "veto", None, scores
        return "accept", proposed, scores

    if (
        total >= min_total
        and best.support >= min_rescue_support
        and best.share >= rescue_share
        and best.anchors >= 1
        and (second <= 0 or best.support >= second * rescue_margin)
    ):
        return "rescue", best.uid, scores
    return "none", None, scores
