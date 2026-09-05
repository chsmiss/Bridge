#!/usr/bin/env python3
"""Stage 16: root-cause experiments for Bridge short-read metagenome assembly.

This stage deliberately avoids metric-target tuning. It isolates two mechanisms
suggested by Stage15 diagnostics:

1. ambiguity-aware approximate read-to-graph threading on the clean k31 graph.
   Sparse k15 minimizer anchors are chained through graph topology; ambiguous
   anchors are allowed when the surrounding anchor chain makes one graph path
   uniquely better. No k-mers or graph edges are added.
2. seed-local rare rescue. All Stage10 cross-k validated rare seeds are used to
   recruit their raw read neighbourhood. Singleton/mercy rescue is enabled only
   inside this targeted pool, never over the full metagenome. Local k17/21/25/31
   assemblies are retained only when they are connected to a trusted seed and
   add raw-supported sequence beyond Stage10.

The outputs are independent plus a combined candidate. Diagnostics focus on
read-path utilization, recruited-library fraction, graph inflation, and evidence
for fresh seed-connected sequence; MetaQUAST/QUAST are downstream observations.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import resource
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import adaptive_k_local_v2 as ak
import graph_path_phaser as gp
import low_abundance_rescue as lr
import repeat_graph_optimizer as rg
import stage10_multik_rescue as s10
import stage14_amplified_methods as s14
import stage15_structural_recovery as s15
import stage789_optimizer as s78


def run(cmd: list[object], *, env: dict[str, str] | None = None) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True, env=env)
    return time.monotonic() - started


@dataclass(frozen=True)
class AnchorEvent:
    pos: int
    candidates: tuple[str, ...]


@dataclass(frozen=True)
class ChainState:
    path: tuple[str, ...]
    last_pos: int
    anchors: int
    ambiguous_anchors: int
    graph_bridges: int
    score: float


@dataclass
class LocalRareEvidence:
    k: int
    name: str
    seq: str
    seed_id: int
    seed_hits: int
    second_seed_hits: int
    fresh31: int
    fresh_fraction: float
    cross_k_sources: int
    cross_k_fraction: float
    raw_supported21: int
    raw_support_fraction: float


class SparseGraphIndex:
    """Oriented k-mer index that drops keys exceeding max occurrences."""

    def __init__(self, graph: gp.Graph, k: int, max_occurrences: int = 12) -> None:
        self.k = k
        members: dict[int, list[str]] = {}
        dropped: set[int] = set()
        for uid, seq in graph.seqs.items():
            seen: set[int] = set()
            for _pos, key in gp.rolling_keys(seq, k):
                if key in seen or key in dropped:
                    continue
                seen.add(key)
                vals = members.setdefault(key, [])
                if uid not in vals:
                    vals.append(uid)
                if len(vals) > max_occurrences:
                    dropped.add(key)
                    members.pop(key, None)
        self.unique: dict[int, str | None] = {}
        self.ambig: dict[int, list[str]] = {}
        for key, vals in members.items():
            if len(vals) == 1:
                self.unique[key] = vals[0]
            elif len(vals) > 1:
                self.unique[key] = None
                self.ambig[key] = vals
        self.dropped_repetitive_keys = len(dropped)


def minimizer_keys(seq: str, k: int, window: int) -> list[tuple[int, int]]:
    kmers = list(gp.rolling_keys(seq, k))
    if not kmers:
        return []
    width = max(1, min(window, len(kmers)))
    result: list[tuple[int, int]] = []
    last: tuple[int, int] | None = None
    for start in range(0, len(kmers) - width + 1):
        chosen = min(kmers[start : start + width], key=lambda item: (item[1], item[0]))
        if chosen != last:
            result.append(chosen)
            last = chosen
    if not result:
        result.append(min(kmers, key=lambda item: (item[1], item[0])))
    return result


def anchor_events(
    seq: str, index: gp.KmerIndex | SparseGraphIndex, window: int
) -> list[AnchorEvent]:
    events: list[AnchorEvent] = []
    for pos, key in minimizer_keys(seq, index.k, window):
        uid = index.unique.get(key)
        if uid is not None:
            events.append(AnchorEvent(pos, (uid,)))
            continue
        cands = index.ambig.get(key)
        if cands:
            events.append(AnchorEvent(pos, tuple(cands)))
    return events


def bridge_transition(
    graph: gp.Graph,
    source: str,
    target: str,
    read_delta: int,
    cache: dict[tuple[str, str], list[str] | None],
    *,
    max_edges: int = 3,
    max_bp: int = 190,
) -> tuple[list[str] | None, bool]:
    if source == target:
        return [source], False
    if target in graph.out.get(source, []):
        return [source, target], False
    key = (source, target)
    if key not in cache:
        cache[key] = s14.unique_bounded_bridge(
            graph,
            source,
            target,
            max_edges=max_edges,
            max_bp=max_bp,
            max_states=240,
        )
    bridge = cache[key]
    if bridge is None:
        return None, False
    internal_bp = sum(
        max(1, len(graph.seqs[uid]) - graph.k) for uid in bridge[1:-1]
    )
    if internal_bp > read_delta + 2 * graph.k + 16:
        return None, False
    return bridge, len(bridge) > 2


def chain_anchor_events(
    events: list[AnchorEvent],
    graph: gp.Graph,
    *,
    beam_width: int = 6,
    cache: dict[tuple[str, str], list[str] | None] | None = None,
) -> ChainState | None:
    if not events:
        return None
    bridge_cache = cache if cache is not None else {}
    beam: list[ChainState] = []
    for event in events:
        candidates = tuple(dict.fromkeys(event.candidates))
        if not candidates:
            continue
        ambiguity = 1 if len(candidates) > 1 else 0
        additions: list[ChainState] = list(beam)
        start_bonus = 1.4 if len(candidates) == 1 else 0.7
        for uid in candidates:
            additions.append(
                ChainState((uid,), event.pos, 1, ambiguity, 0, start_bonus)
            )
        for state in beam:
            if event.pos <= state.last_pos:
                continue
            delta = event.pos - state.last_pos
            last = state.path[-1]
            for uid in candidates:
                bridge, bridged = bridge_transition(
                    graph, last, uid, delta, bridge_cache
                )
                if bridge is None:
                    continue
                if uid == last:
                    path = state.path
                    transition_bonus = 0.15
                else:
                    extra = tuple(bridge[1:])
                    if any(
                        node in state.path
                        or graph.rev.get(node, node) in state.path
                        for node in extra
                    ):
                        continue
                    path = state.path + extra
                    transition_bonus = 0.9 if not bridged else 0.55
                anchor_bonus = 1.5 if len(candidates) == 1 else 0.9
                additions.append(
                    ChainState(
                        path,
                        event.pos,
                        state.anchors + 1,
                        state.ambiguous_anchors + ambiguity,
                        state.graph_bridges + int(bridged),
                        state.score + anchor_bonus + transition_bonus,
                    )
                )
        best_by_signature: dict[tuple[tuple[str, ...], int], ChainState] = {}
        for state in additions:
            sig = (state.path, state.last_pos)
            old = best_by_signature.get(sig)
            if old is None or (state.anchors, state.score) > (
                old.anchors,
                old.score,
            ):
                best_by_signature[sig] = state
        beam = sorted(
            best_by_signature.values(),
            key=lambda state: (
                -state.anchors,
                -state.score,
                -len(state.path),
                state.path,
            ),
        )[:beam_width]

    eligible = [
        state for state in beam if state.anchors >= 2 and len(state.path) >= 2
    ]
    if not eligible:
        return None
    eligible.sort(
        key=lambda state: (
            -state.anchors,
            -state.score,
            -len(state.path),
            state.path,
        )
    )
    best = eligible[0]
    second = next(
        (state for state in eligible[1:] if state.path != best.path), None
    )
    if second is not None:
        if best.anchors < second.anchors:
            return None
        if best.anchors == second.anchors and best.score < second.score + 0.9:
            return None
    return best


def chain_sequence(
    seq: str,
    graph: gp.Graph,
    index: gp.KmerIndex | SparseGraphIndex,
    *,
    window: int = 8,
    beam_width: int = 6,
    cache: dict[tuple[str, str], list[str] | None] | None = None,
) -> tuple[ChainState | None, int]:
    events = anchor_events(seq, index, window)
    return chain_anchor_events(
        events, graph, beam_width=beam_width, cache=cache
    ), sum(len(event.candidates) > 1 for event in events)


def collect_chained_contexts(
    graph: gp.Graph,
    index: gp.KmerIndex | SparseGraphIndex,
    read1: Path,
    read2: Path,
    *,
    exact_index: gp.KmerIndex | None = None,
    max_context: int = 10,
) -> tuple[Counter[tuple[str, ...]], dict[str, int]]:
    contexts: Counter[tuple[str, ...]] = Counter()
    stats = {
        "reads": 0,
        "reads_with_indexed_anchors": 0,
        "reads_with_ambiguous_anchors": 0,
        "threaded_reads": 0,
        "approx_only_threaded_reads": 0,
        "threaded_with_ambiguous_anchor": 0,
        "graph_bridge_reads": 0,
        "anchors_on_accepted_chains": 0,
        "ambiguous_anchors_on_accepted_chains": 0,
        "contexts": 0,
    }
    cache: dict[tuple[str, str], list[str] | None] = {}
    for fastq in (read1, read2):
        for _name, seq in gp.read_fastq(fastq):
            stats["reads"] += 1
            events = anchor_events(seq, index, 8)
            if events:
                stats["reads_with_indexed_anchors"] += 1
            ambiguous_events = sum(
                len(event.candidates) > 1 for event in events
            )
            if ambiguous_events:
                stats["reads_with_ambiguous_anchors"] += 1
            chain = chain_anchor_events(
                events, graph, beam_width=6, cache=cache
            )
            if chain is None:
                continue
            stats["threaded_reads"] += 1
            if exact_index is not None and not gp.thread_sequence(
                seq, graph, exact_index, None
            ):
                stats["approx_only_threaded_reads"] += 1
            stats["anchors_on_accepted_chains"] += chain.anchors
            stats[
                "ambiguous_anchors_on_accepted_chains"
            ] += chain.ambiguous_anchors
            if chain.ambiguous_anchors:
                stats["threaded_with_ambiguous_anchor"] += 1
            if chain.graph_bridges:
                stats["graph_bridge_reads"] += 1
            gp.add_context(contexts, list(chain.path), max_context)
    stats["contexts"] = len(contexts)
    stats["bridge_cache_entries"] = len(cache)
    return contexts, stats


def conservative_increment(
    baseline: Counter[tuple[str, ...]],
    approximate: Counter[tuple[str, ...]],
    *,
    max_weight: int = 6,
) -> Counter[tuple[str, ...]]:
    result: Counter[tuple[str, ...]] = Counter()
    for key, support in approximate.items():
        delta = max(0, support - baseline.get(key, 0))
        minimum = 1 if len(key) >= 3 else 2
        if delta >= minimum:
            result[key] = min(max_weight, delta)
    return result


def build_approx_thread_candidate(
    scripts: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    read1: Path,
    read2: Path,
    timings: dict[str, float],
) -> tuple[Path, dict[str, object]]:
    base = pipeline_dir / "current_pipeline"
    graph_opt = pipeline_dir / "graph_optimizer"
    repeat_opt = pipeline_dir / "repeat_optimizer"
    target_gfa = base / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = (
        base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    )
    projection_haplotigs = (
        base / "iterative" / "k21_recall" / "haplotigs.fasta"
    )
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"
    base_paths = graph_opt / "stage4_second_pass.paths.tsv"

    graph = gp.Graph.from_gfa(target_gfa)
    membership = gp.preliminary_membership(rg.load_paths(base_paths))
    exact_index = gp.KmerIndex(graph, 31)
    sparse_index = SparseGraphIndex(graph, 15, max_occurrences=6)

    started = time.monotonic()
    exact_ctx, exact_stats = gp.collect_read_contexts(
        graph, exact_index, read1, read2, None, 10
    )
    approx_ctx, approx_stats = collect_chained_contexts(
        graph,
        sparse_index,
        read1,
        read2,
        exact_index=exact_index,
        max_context=10,
    )
    increment = conservative_increment(exact_ctx, approx_ctx)
    enhanced_raw = Counter(exact_ctx)
    enhanced_raw.update(increment)
    timings["approx_thread_full_library"] = time.monotonic() - started

    proj_ctx, high_ctx, projection_stats = rg.collect_projection_contexts(
        graph,
        exact_index,
        [projection_primary, projection_haplotigs],
        [highk_gfa],
        repeat_opt,
        8,
    )
    second_ctx, second_stats = gp.collect_read_contexts(
        graph, exact_index, read1, read2, membership, 8
    )
    for key in list(second_ctx):
        base_support = exact_ctx.get(key, 0)
        if second_ctx[key] <= base_support:
            del second_ctx[key]
        else:
            second_ctx[key] -= base_support
    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph, exact_index, read1, read2, membership, 8, 8, 420
    )
    repeat_ctx = rg.combined_contexts(second_ctx, pair_ctx)
    all_ctx = rg.combined_contexts(
        enhanced_raw, proj_ctx, high_ctx, repeat_ctx
    )
    simplified, simplify_stats = rg.simplify_graph(graph, all_ctx)
    paths, resolve_stats = s78.resolve_lookahead_seeded_paths(
        simplified,
        enhanced_raw,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        0.70,
        4,
        200,
        4,
        5,
        0.70,
        0.58,
        1.10,
    )

    outdir = pipeline_dir / "stage16_root_cause" / "approx_thread"
    outdir.mkdir(parents=True, exist_ok=True)
    raw_fasta = outdir / "approx_thread_raw.fasta"
    write_stats = gp.write_paths(
        paths,
        simplified,
        raw_fasta,
        outdir / "approx_thread.paths.tsv",
        200,
    )
    final = s78.emit_stage(
        scripts,
        raw_fasta,
        strict_baseline,
        pipeline_dir,
        "stage16_approx_thread",
        31,
        timings,
    )
    copied = (
        outdir / "candidate_approx_thread" / "primary_contigs.fasta"
    )
    copied.parent.mkdir(parents=True, exist_ok=True)
    copied.write_bytes(final.read_bytes())
    return copied, {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "exact_k31": exact_stats,
        "approx_k15_minimizer": {
            **approx_stats,
            "dropped_repetitive_keys": sparse_index.dropped_repetitive_keys,
        },
        "incremental_contexts": len(increment),
        "incremental_context_weight": sum(increment.values()),
        "projection": projection_stats,
        "second_pass": second_stats,
        "pair_threading": pair_stats,
        "simplification": simplify_stats,
        "path_resolution": resolve_stats,
        "write": write_stats,
    }


def build_all_seed_signatures(
    seed_fasta: Path,
    baseline: Path,
    *,
    k: int = 21,
    min_signature_kmers: int = 4,
) -> tuple[
    list[tuple[str, str]], dict[str, int], dict[int, set[str]], dict[str, int]
]:
    baseline_kmers, _ = lr.backbone_kmers(baseline, k)
    raw_seeds = [
        (name, lr.canonical(seq))
        for name, seq in lr.fasta_records(seed_fasta)
    ]
    memberships: dict[str, set[int]] = defaultdict(set)
    for sid, (_name, seq) in enumerate(raw_seeds):
        novel = set(lr.kmers(seq, k)) - baseline_kmers
        for mer in novel:
            memberships[mer].add(sid)
    unique_by_old: dict[int, set[str]] = defaultdict(set)
    for mer, sids in memberships.items():
        if len(sids) == 1:
            unique_by_old[next(iter(sids))].add(mer)

    kept_old = [
        sid
        for sid in range(len(raw_seeds))
        if len(unique_by_old[sid]) >= min_signature_kmers
    ]
    old_to_new = {old: new for new, old in enumerate(kept_old)}
    seeds = [raw_seeds[old] for old in kept_old]
    signatures: dict[str, int] = {}
    signature_sets: dict[int, set[str]] = {}
    for old, new in old_to_new.items():
        signature_sets[new] = unique_by_old[old]
        for mer in unique_by_old[old]:
            signatures[mer] = new
    return seeds, signatures, signature_sets, {
        "input_seed_records": len(raw_seeds),
        "signature_eligible_seeds": len(seeds),
        "initial_unique_signatures": len(signatures),
        "seeds_without_enough_unique_signature": len(raw_seeds) - len(seeds),
    }


def score_pair(
    seq1: str,
    seq2: str,
    signatures: dict[str, int],
    k: int,
    stride: int,
) -> Counter[int]:
    scores: Counter[int] = Counter()
    for seq in (seq1, seq2):
        for mer in lr.kmers(seq, k, stride):
            sid = signatures.get(mer)
            if sid is not None:
                scores[sid] += 1
    return scores


def unique_assignment(
    scores: Counter[int], min_hits: int, margin: int
) -> tuple[int | None, bool]:
    ranked = scores.most_common(2)
    if not ranked or ranked[0][1] < min_hits:
        return None, False
    ambiguous = (
        len(ranked) > 1 and ranked[0][1] - ranked[1][1] < margin
    )
    if ambiguous:
        return None, True
    return ranked[0][0], False


def recruit_all_seed_pairs(
    read1: Path,
    read2: Path,
    seeds: list[tuple[str, str]],
    initial_signatures: dict[str, int],
    baseline: Path,
) -> tuple[dict[int, int], dict[str, int], dict[str, int]]:
    pair_to_seed: dict[int, int] = {}
    initial_rejected_ambiguous = 0
    total_pairs = 0
    for idx, (left, right) in enumerate(
        zip(ak.fastq_records(read1), ak.fastq_records(read2))
    ):
        total_pairs += 1
        sid, ambiguous = unique_assignment(
            score_pair(
                left[1], right[1], initial_signatures, 21, 2
            ),
            2,
            1,
        )
        if ambiguous:
            initial_rejected_ambiguous += 1
        if sid is not None:
            pair_to_seed[idx] = sid
    initial_pairs = len(pair_to_seed)

    baseline17, _ = lr.backbone_kmers(baseline, 17)
    counts_by_seed: dict[int, Counter[str]] = defaultdict(Counter)
    for idx, (left, right) in enumerate(
        zip(ak.fastq_records(read1), ak.fastq_records(read2))
    ):
        sid = pair_to_seed.get(idx)
        if sid is None:
            continue
        fragment_mers = set(lr.kmers(left[1], 17, 2)) | set(
            lr.kmers(right[1], 17, 2)
        )
        for mer in fragment_mers:
            if mer not in baseline17:
                counts_by_seed[sid][mer] += 1

    mer_membership: dict[str, set[int]] = defaultdict(set)
    for sid, counts in counts_by_seed.items():
        for mer, count in counts.items():
            if count >= 2:
                mer_membership[mer].add(sid)
    expanded = {
        mer: next(iter(sids))
        for mer, sids in mer_membership.items()
        if len(sids) == 1
    }

    second_rejected_ambiguous = 0
    second_pairs = 0
    for idx, (left, right) in enumerate(
        zip(ak.fastq_records(read1), ak.fastq_records(read2))
    ):
        if idx in pair_to_seed:
            continue
        sid, ambiguous = unique_assignment(
            score_pair(left[1], right[1], expanded, 17, 2), 3, 2
        )
        if ambiguous:
            second_rejected_ambiguous += 1
        if sid is not None:
            pair_to_seed[idx] = sid
            second_pairs += 1

    per_seed = Counter(pair_to_seed.values())
    stats = {
        "total_library_pairs": total_pairs,
        "seed_loci": len(seeds),
        "initial_pairs": initial_pairs,
        "second_round_pairs": second_pairs,
        "total_recruited_pairs": len(pair_to_seed),
        "expanded_unique_signatures": len(expanded),
        "initial_ambiguous_pair_rejects": initial_rejected_ambiguous,
        "second_ambiguous_pair_rejects": second_rejected_ambiguous,
        "seeds_with_recruited_pairs": len(per_seed),
        "max_pairs_per_seed": max(per_seed.values(), default=0),
        "median_pairs_per_seed": (
            int(sorted(per_seed.values())[len(per_seed) // 2])
            if per_seed
            else 0
        ),
    }
    return pair_to_seed, stats, expanded


def write_targeted_pool(
    read1: Path,
    read2: Path,
    pair_to_seed: dict[int, int],
    out1: Path,
    out2: Path,
) -> None:
    out1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out1, "wt") as left_out, gzip.open(
        out2, "wt"
    ) as right_out:
        for idx, (left, right) in enumerate(
            zip(ak.fastq_records(read1), ak.fastq_records(read2))
        ):
            if idx not in pair_to_seed:
                continue
            ak.write_record(left_out, left)
            ak.write_record(right_out, right)


def targeted_raw_kmer_support(
    read1: Path, read2: Path, k: int = 21
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for left, right in zip(
        ak.fastq_records(read1), ak.fastq_records(read2)
    ):
        fragment = set(lr.kmers(left[1], k, 1)) | set(
            lr.kmers(right[1], k, 1)
        )
        counts.update(fragment)
    return counts


def assemble_targeted_k(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    outdir: Path,
    k: int,
    threads: int,
) -> float:
    mercy = {17: 72, 21: 56, 25: 44, 31: 32}.get(k, 32)
    env = os.environ.copy()
    env["BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION"] = "0.80"
    env["BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY"] = "35"
    env["BRIDGEASM_MATE_TERMINAL_MERCY_KMERS"] = "96"
    return run(
        [
            bridgeasm,
            "assemble",
            "-1",
            read1,
            "-2",
            read2,
            "-o",
            outdir,
            "-k",
            k,
            "--min-count",
            2,
            "--mercy-max-kmers",
            mercy,
            "--mercy-min-support",
            1,
            "--mercy-min-quality",
            30,
            "--min-read-support",
            2,
            "--min-pair-support",
            2,
            "--min-primary-support",
            3,
            "--primary-dominance",
            0.82,
            "--threaded-path-cover",
            "--major-path-cover",
            "--path-cover-secondary-dominance",
            0.20,
            "--min-contig-length",
            200,
            "--threads",
            threads,
        ],
        env=env,
    )


def local_evidence(
    assembly_inputs: dict[int, list[Path]],
    seeds: list[tuple[str, str]],
    strict_baseline: Path,
    raw21_counts: Counter[str],
) -> list[LocalRareEvidence]:
    baseline31, _ = lr.backbone_kmers(strict_baseline, 31)
    baseline21, _ = lr.backbone_kmers(strict_baseline, 21)
    seed_membership: dict[str, set[int]] = defaultdict(set)
    for sid, (_name, seq) in enumerate(seeds):
        for mer in set(lr.kmers(seq, 21)):
            seed_membership[mer].add(sid)

    pools: dict[int, set[str]] = defaultdict(set)
    raw_records: list[tuple[int, str, str]] = []
    for k, paths in sorted(assembly_inputs.items()):
        for path in paths:
            if not path.exists():
                continue
            for name, seq0 in lr.fasta_records(path):
                seq = lr.canonical(seq0)
                if len(seq) < 200:
                    continue
                raw_records.append((k, name, seq))
                pools[k].update(lr.kmers(seq, 21))

    result: list[LocalRareEvidence] = []
    for k, name, seq in raw_records:
        all31 = set(lr.kmers(seq, 31))
        fresh31 = all31 - baseline31
        if not fresh31:
            continue
        all21 = set(lr.kmers(seq, 21))
        seed_scores: Counter[int] = Counter()
        for mer in all21:
            for sid in seed_membership.get(mer, ()):
                seed_scores[sid] += 1
        ranked = seed_scores.most_common(2)
        if not ranked:
            continue
        sid, seed_hits = ranked[0]
        second_hits = ranked[1][1] if len(ranked) > 1 else 0
        if seed_hits < 4 or seed_hits < second_hits + 2:
            continue

        fresh21 = all21 - baseline21
        cross_sources = 0
        union_other: set[str] = set()
        for other_k, pool in pools.items():
            if other_k == k:
                continue
            union_other.update(pool)
            if fresh21 and len(fresh21 & pool) / len(fresh21) >= 0.20:
                cross_sources += 1
        cross_fraction = (
            len(fresh21 & union_other) / len(fresh21)
            if fresh21
            else 0.0
        )
        raw_supported = sum(
            raw21_counts.get(mer, 0) >= 2 for mer in fresh21
        )
        raw_fraction = (
            raw_supported / len(fresh21) if fresh21 else 0.0
        )
        result.append(
            LocalRareEvidence(
                k=k,
                name=name,
                seq=seq,
                seed_id=sid,
                seed_hits=seed_hits,
                second_seed_hits=second_hits,
                fresh31=len(fresh31),
                fresh_fraction=len(fresh31) / max(1, len(all31)),
                cross_k_sources=cross_sources,
                cross_k_fraction=cross_fraction,
                raw_supported21=raw_supported,
                raw_support_fraction=raw_fraction,
            )
        )
    return result


def select_local_evidence(
    items: list[LocalRareEvidence], strict_baseline: Path
) -> list[LocalRareEvidence]:
    baseline31, baseline_bases = lr.backbone_kmers(
        strict_baseline, 31
    )
    items = sorted(
        items,
        key=lambda item: (
            item.seed_id,
            -item.cross_k_sources,
            -item.cross_k_fraction,
            -item.raw_support_fraction,
            -item.raw_supported21,
            -item.fresh31,
            -len(item.seq),
            item.k,
            item.seq,
        ),
    )
    selected: list[LocalRareEvidence] = []
    selected_fresh: set[str] = set()
    per_seed: Counter[int] = Counter()
    total_bases = 0
    max_bases = max(100_000, int(baseline_bases * 0.12))
    for item in items:
        evidence_ok = (
            item.cross_k_sources >= 1
            or item.cross_k_fraction >= 0.45
            or (
                item.raw_supported21 >= 8
                and item.raw_support_fraction >= 0.45
            )
        )
        if not evidence_ok:
            continue
        novel = set(lr.kmers(item.seq, 31)) - baseline31
        fresh = novel - selected_fresh
        if len(fresh) < max(8, math.ceil(0.25 * len(novel))):
            continue
        if per_seed[item.seed_id] >= 3:
            continue
        if total_bases + len(item.seq) > max_bases:
            continue
        selected.append(item)
        selected_fresh.update(fresh)
        per_seed[item.seed_id] += 1
        total_bases += len(item.seq)
    return selected


def write_local_evidence(
    items: Iterable[LocalRareEvidence], path: Path
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "k\tname\tlength\tseed_id\tseed_hits\tsecond_seed_hits"
            "\tfresh31\tfresh_fraction\tcross_k_sources\tcross_k_fraction"
            "\traw_supported21\traw_support_fraction\n"
        )
        for item in items:
            handle.write(
                f"{item.k}\t{item.name}\t{len(item.seq)}"
                f"\t{item.seed_id}\t{item.seed_hits}"
                f"\t{item.second_seed_hits}\t{item.fresh31}"
                f"\t{item.fresh_fraction:.6f}"
                f"\t{item.cross_k_sources}"
                f"\t{item.cross_k_fraction:.6f}"
                f"\t{item.raw_supported21}"
                f"\t{item.raw_support_fraction:.6f}\n"
            )


def profile(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.exists() else {}


def build_local_rare_candidate(
    scripts: Path,
    bridgeasm: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    backbone: Path,
    read1: Path,
    read2: Path,
    threads: int,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    stage10 = pipeline_dir / "stage10_multik_rescue"
    seed_fasta = stage10 / "multik_strict_additions.fasta"
    outdir = pipeline_dir / "stage16_root_cause" / "local_rare"
    outdir.mkdir(parents=True, exist_ok=True)

    seeds, signatures, _signature_sets, signature_stats = (
        build_all_seed_signatures(
            seed_fasta, backbone, k=21, min_signature_kmers=4
        )
    )
    started = time.monotonic()
    pair_to_seed, recruit_stats, expanded = recruit_all_seed_pairs(
        read1, read2, seeds, signatures, backbone
    )
    raw1 = outdir / "targeted_raw_R1.fastq.gz"
    raw2 = outdir / "targeted_raw_R2.fastq.gz"
    pair_counts = Counter(pair_to_seed.values())
    with (outdir / "seed_recruitment.tsv").open("w") as handle:
        handle.write("seed_id\tseed_name\trecruited_pairs\n")
        for sid, (name, _seq) in enumerate(seeds):
            handle.write(
                f"{sid}\t{name}\t{pair_counts.get(sid, 0)}\n"
            )
    write_targeted_pool(read1, read2, pair_to_seed, raw1, raw2)
    timings["local_rare_recruitment"] = time.monotonic() - started

    seed_v1 = outdir / "trusted_seed_R1.fastq.gz"
    seed_v2 = outdir / "trusted_seed_R2.fastq.gz"
    timings["local_rare_seed_virtualization"] = s15.virtualize(
        scripts,
        [seed_fasta],
        seed_v1,
        seed_v2,
        read_length=91,
        insert_size=190,
        stride=60,
        min_length=190,
    )
    aug1 = outdir / "targeted_aug_R1.fastq.gz"
    aug2 = outdir / "targeted_aug_R2.fastq.gz"
    s15.concat_gzip([raw1, seed_v1], aug1)
    s15.concat_gzip([raw2, seed_v2], aug2)
    raw21_counts = targeted_raw_kmer_support(raw1, raw2, 21)

    assembly_inputs: dict[int, list[Path]] = {}
    profiles: dict[str, object] = {}
    for k in (17, 21, 25, 31):
        asm = outdir / f"k{k}"
        timings[f"local_rare_k{k}"] = assemble_targeted_k(
            bridgeasm, aug1, aug2, asm, k, threads
        )
        assembly_inputs[k] = [
            asm / "primary_contigs.fasta",
            asm / "haplotigs.fasta",
        ]
        profiles[f"k{k}"] = profile(asm / "run_profile.json")

    evidence = local_evidence(
        assembly_inputs, seeds, strict_baseline, raw21_counts
    )
    selected = select_local_evidence(evidence, strict_baseline)
    write_local_evidence(evidence, outdir / "local_rare_evidence.tsv")
    write_local_evidence(
        selected, outdir / "local_rare_selected.tsv"
    )
    additions = outdir / "local_rare_additions.fasta"
    s14.write_fasta(
        (
            (
                f"stage16_seed{item.seed_id}_k{item.k}_{idx:06d}",
                item.seq,
            )
            for idx, item in enumerate(selected, 1)
        ),
        additions,
    )
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_local_rare",
        timings,
        min_overlap=81,
    )
    selected_seeds = {item.seed_id for item in selected}
    return final, additions, {
        **signature_stats,
        **recruit_stats,
        "recruited_library_fraction": len(pair_to_seed)
        / max(1, recruit_stats["total_library_pairs"]),
        "expanded_signatures": len(expanded),
        "raw_supported_k21": len(raw21_counts),
        "evidence_records": len(evidence),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "selected_seed_loci": len(selected_seeds),
        "selected_fresh31": sum(item.fresh31 for item in selected),
        "assembly_profiles": profiles,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    started = time.monotonic()
    scripts = Path(__file__).resolve().parent
    pipeline = args.pipeline_dir
    stage10 = pipeline / "stage10_multik_rescue"
    strict_baseline = (
        stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    )
    backbone = pipeline / "bridge_backbone.fasta"
    required = [
        args.bridgeasm,
        args.read1,
        args.read2,
        strict_baseline,
        backbone,
        stage10 / "multik_strict_additions.fasta",
        pipeline
        / "current_pipeline"
        / "iterative"
        / "k31_resolve"
        / "assembly.gfa",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit(
            "missing Stage16 inputs: " + ", ".join(missing)
        )

    root = pipeline / "stage16_root_cause"
    root.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    approx_final, approx_stats = build_approx_thread_candidate(
        scripts,
        pipeline,
        strict_baseline,
        args.read1,
        args.read2,
        timings,
    )
    local_final, local_additions, local_stats = (
        build_local_rare_candidate(
            scripts,
            args.bridgeasm,
            pipeline,
            strict_baseline,
            backbone,
            args.read1,
            args.read2,
            args.threads,
            timings,
        )
    )
    combined = s14.make_bridge_candidate(
        scripts,
        approx_final,
        [local_additions],
        root / "combined" / "candidate_combined",
        timings,
        min_overlap=81,
    )

    stats = {
        "pipeline": "bridge-stage16-root-cause-v1",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "metric_targets": False,
            "approx_thread": (
                "k15 sparse minimizer anchors + ambiguity-aware graph "
                "chaining on clean k31 graph"
            ),
            "local_rare": (
                "all signature-eligible Stage10 seeds recruit raw pairs; "
                "sparse rescue only in targeted pool"
            ),
            "local_validation": (
                "trusted-seed connection + fresh Stage10-novel sequence + "
                "raw/cross-k evidence"
            ),
            "sequence_join": "exact overlap >=81 bp",
        },
        "methods": {
            "approx_thread": approx_stats,
            "local_rare": local_stats,
        },
        "outputs": {
            "approx_thread": str(approx_final),
            "local_rare": str(local_final),
            "combined": str(combined),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss,
    }
    (root / "stage16_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
