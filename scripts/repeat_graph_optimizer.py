#!/usr/bin/env python3
"""Conservative repeat traversal and primary-graph simplification.

This module is intentionally downstream of full-read phasing. It adds two
cumulative candidate backbones:
  5. mate-pair repeat traversal: combine graph-threaded mates into longer
     contexts and re-resolve ambiguous exits.
  6. conservative graph simplification: mask only weak high-confidence tips or
     dominated branch edges that lack long-context support, then re-resolve.

No new sequence adjacency is invented; every emitted transition exists in the
input target GFA.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import graph_path_phaser as gp


def load_paths(path: Path) -> list[list[str]]:
    paths: list[list[str]] = []
    with path.open() as handle:
        header = next(handle, None)
        if header is None:
            return paths
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            nodes = [node for node in fields[3].split(",") if node]
            if nodes:
                paths.append(nodes)
    return paths


def reverse_path(path: list[str], graph: gp.Graph) -> list[str]:
    return [graph.rev.get(uid, uid) for uid in reversed(path)]


def best_segment(segments: list[list[str]], graph: gp.Graph) -> list[str] | None:
    if not segments:
        return None
    return max(
        segments,
        key=lambda path: (
            len(path),
            sum(max(1, len(graph.seqs[uid]) - graph.k) for uid in path),
        ),
    )


def join_paths(
    left: list[str],
    right: list[str],
    graph: gp.Graph,
    max_edges: int,
    max_span: int,
) -> list[str] | None:
    if not left or not right:
        return None
    if left[-1] == right[0]:
        merged = left + right[1:]
        return merged if len(set(merged)) == len(merged) else None
    bridge = gp.unique_short_bridge(graph, left[-1], right[0], max_edges)
    if bridge is None:
        return None
    bridge_bp = sum(max(1, len(graph.seqs[uid]) - graph.k) for uid in bridge[1:])
    if bridge_bp > max_span:
        return None
    merged = left + bridge[1:-1] + right
    return merged if len(set(merged)) == len(merged) else None


def collect_pair_contexts(
    graph: gp.Graph,
    index: gp.KmerIndex,
    read1: Path,
    read2: Path,
    membership: dict[str, list[tuple[int, int]]],
    max_context: int,
    max_bridge_edges: int,
    max_pair_span: int,
) -> tuple[Counter[tuple[str, ...]], dict[str, int]]:
    contexts: Counter[tuple[str, ...]] = Counter()
    stats = {
        "pairs": 0,
        "both_threaded": 0,
        "joined_pairs": 0,
        "contexts": 0,
        "support": 0,
    }
    it1 = gp.read_fastq(read1)
    it2 = gp.read_fastq(read2)
    for left_rec, right_rec in zip(it1, it2):
        stats["pairs"] += 1
        left = best_segment(
            gp.thread_sequence(left_rec[1], graph, index, membership), graph
        )
        right = best_segment(
            gp.thread_sequence(right_rec[1], graph, index, membership), graph
        )
        if left is None or right is None:
            continue
        stats["both_threaded"] += 1
        candidates: list[list[str]] = []
        for a, b in (
            (left, reverse_path(right, graph)),
            (right, reverse_path(left, graph)),
            (left, right),
            (right, left),
        ):
            merged = join_paths(a, b, graph, max_bridge_edges, max_pair_span)
            if merged is not None:
                candidates.append(merged)
        if not candidates:
            continue
        candidates.sort(key=lambda path: (-len(path), tuple(path)))
        best = candidates[0]
        if (
            len(candidates) > 1
            and candidates[1] != best
            and len(candidates[1]) == len(best)
        ):
            continue
        stats["joined_pairs"] += 1
        gp.add_context(contexts, best, max_context, weight=3)
    stats["contexts"] = len(contexts)
    stats["support"] = sum(contexts.values())
    return contexts, stats


def collect_projection_contexts(
    graph: gp.Graph,
    index: gp.KmerIndex,
    projections: list[Path],
    highk_gfas: list[Path],
    out: Path,
    max_context: int,
) -> tuple[Counter[tuple[str, ...]], Counter[tuple[str, ...]], dict[str, object]]:
    proj_ctx, proj_stats = gp.collect_fasta_contexts(
        graph, index, projections, max_context
    )
    high_ctx: Counter[tuple[str, ...]] = Counter()
    high_stats = {
        "records": 0,
        "segments": 0,
        "accepted_segments": 0,
        "contexts": 0,
    }
    for high_gfa in highk_gfas:
        high_graph = gp.Graph.from_gfa(high_gfa)
        tmp = out / (high_gfa.parent.name + ".repeat_highk_unitigs.fasta")
        with tmp.open("w") as handle:
            for uid, seq in high_graph.seqs.items():
                handle.write(f">{uid}\n{seq}\n")
        ctx, stage_stats = gp.collect_fasta_contexts(
            graph, index, [tmp], max_context, only_ambiguous=True
        )
        high_ctx.update(ctx)
        for key in ("records", "segments", "accepted_segments"):
            high_stats[key] += stage_stats.get(key, 0)
    high_stats["contexts"] = len(high_ctx)
    return proj_ctx, high_ctx, {"projection": proj_stats, "highk": high_stats}


def combined_contexts(*counters: Counter[tuple[str, ...]]) -> Counter[tuple[str, ...]]:
    result: Counter[tuple[str, ...]] = Counter()
    for counter in counters:
        result.update(counter)
    return result


def long_edge_support(counter: Counter[tuple[str, ...]]) -> Counter[tuple[str, str]]:
    support: Counter[tuple[str, str]] = Counter()
    for context, count in counter.items():
        if len(context) < 3:
            continue
        for edge in zip(context, context[1:]):
            support[edge] += count
    return support


def edge_score(
    graph: gp.Graph,
    edge: tuple[str, str],
    context: Counter[tuple[str, ...]],
    long_support: Counter[tuple[str, str]],
) -> float:
    ev = graph.edge.get(edge, gp.EdgeEvidence())
    return (
        ev.direct * 2.0
        + ev.gapped
        + ev.pairs * 4.0
        + context.get(edge, 0) * 2.0
        + long_support.get(edge, 0) * 3.0
    )


def copy_without_edges(graph: gp.Graph, masked: set[tuple[str, str]]) -> gp.Graph:
    result = gp.Graph()
    result.k = graph.k
    result.seqs = dict(graph.seqs)
    result.coverage = dict(graph.coverage)
    result.rev = dict(graph.rev)
    for (src, dst), evidence in graph.edge.items():
        if (src, dst) in masked:
            continue
        result.edge[(src, dst)] = evidence
        result.out[src].append(dst)
        result.inc[dst].append(src)
    for uid in result.seqs:
        result.out[uid] = sorted(set(result.out.get(uid, [])))
        result.inc[uid] = sorted(set(result.inc.get(uid, [])))
    return result


def simplify_graph(
    graph: gp.Graph,
    context: Counter[tuple[str, ...]],
    max_tip_bp: int = 400,
    dominance: float = 0.82,
    max_loser_fraction: float = 0.22,
) -> tuple[gp.Graph, dict[str, int]]:
    long_support = long_edge_support(context)
    masked: set[tuple[str, str]] = set()
    tip_edges = 0
    dominated_edges = 0

    for src, children in list(graph.out.items()):
        if len(children) <= 1:
            continue
        src_cov = max(1e-6, graph.coverage.get(src, 0.0))
        for dst in children:
            ev = graph.edge.get((src, dst), gp.EdgeEvidence())
            dst_cov = graph.coverage.get(dst, 0.0)
            if (
                len(graph.out.get(dst, [])) == 0
                and len(graph.seqs[dst]) <= max_tip_bp
                and ev.direct <= 1
                and ev.pairs == 0
                and long_support.get((src, dst), 0) == 0
                and dst_cov <= 0.35 * src_cov
            ):
                masked.add((src, dst))
                tip_edges += 1

    for dst, parents in list(graph.inc.items()):
        if len(parents) <= 1:
            continue
        dst_cov = max(1e-6, graph.coverage.get(dst, 0.0))
        for src in parents:
            ev = graph.edge.get((src, dst), gp.EdgeEvidence())
            src_cov = graph.coverage.get(src, 0.0)
            if (
                len(graph.inc.get(src, [])) == 0
                and len(graph.seqs[src]) <= max_tip_bp
                and ev.direct <= 1
                and ev.pairs == 0
                and long_support.get((src, dst), 0) == 0
                and src_cov <= 0.35 * dst_cov
            ):
                masked.add((src, dst))
                tip_edges += 1

    for src, children in list(graph.out.items()):
        active = [dst for dst in children if (src, dst) not in masked]
        if len(active) <= 1:
            continue
        ranked = sorted(
            (
                (edge_score(graph, (src, dst), context, long_support), dst)
                for dst in active
            ),
            reverse=True,
        )
        best_score, best_dst = ranked[0]
        total = sum(max(0.0, score) for score, _ in ranked)
        if best_score <= 0 or total <= 0 or best_score / total < dominance:
            continue
        best_long = long_support.get((src, best_dst), 0)
        best_ev = graph.edge.get((src, best_dst), gp.EdgeEvidence())
        if best_long < 2 and best_ev.direct < 6 and best_ev.pairs < 2:
            continue
        best_cov = max(1e-6, graph.coverage.get(best_dst, 0.0))
        for score, dst in ranked[1:]:
            if score > best_score * max_loser_fraction:
                continue
            ev = graph.edge.get((src, dst), gp.EdgeEvidence())
            if long_support.get((src, dst), 0) > 0:
                continue
            if ev.pairs > 0 or ev.direct > max(2, best_ev.direct // 3):
                continue
            dst_cov = graph.coverage.get(dst, 0.0)
            if dst_cov > 0.55 * best_cov:
                continue
            masked.add((src, dst))
            dominated_edges += 1

    for dst, parents in list(graph.inc.items()):
        active = [src for src in parents if (src, dst) not in masked]
        if len(active) <= 1:
            continue
        ranked = sorted(
            (
                (edge_score(graph, (src, dst), context, long_support), src)
                for src in active
            ),
            reverse=True,
        )
        best_score, best_src = ranked[0]
        total = sum(max(0.0, score) for score, _ in ranked)
        if best_score <= 0 or total <= 0 or best_score / total < dominance:
            continue
        best_long = long_support.get((best_src, dst), 0)
        best_ev = graph.edge.get((best_src, dst), gp.EdgeEvidence())
        if best_long < 2 and best_ev.direct < 6 and best_ev.pairs < 2:
            continue
        best_cov = max(1e-6, graph.coverage.get(best_src, 0.0))
        for score, src in ranked[1:]:
            if score > best_score * max_loser_fraction:
                continue
            ev = graph.edge.get((src, dst), gp.EdgeEvidence())
            if long_support.get((src, dst), 0) > 0:
                continue
            if ev.pairs > 0 or ev.direct > max(2, best_ev.direct // 3):
                continue
            src_cov = graph.coverage.get(src, 0.0)
            if src_cov > 0.55 * best_cov:
                continue
            masked.add((src, dst))
            dominated_edges += 1

    simplified = copy_without_edges(graph, masked)
    return simplified, {
        "masked_edges": len(masked),
        "tip_edges": tip_edges,
        "dominated_edges": dominated_edges,
        "ambiguous_nodes_before": sum(graph.ambiguous(uid) for uid in graph.seqs),
        "ambiguous_nodes_after": sum(
            simplified.ambiguous(uid) for uid in simplified.seqs
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--base-paths", type=Path, required=True)
    ap.add_argument("--projection", type=Path, action="append", default=[])
    ap.add_argument("--highk-gfa", type=Path, action="append", default=[])
    ap.add_argument("-o", "--output-dir", type=Path, required=True)
    ap.add_argument("--anchor-k", type=int, default=31)
    ap.add_argument("--max-context", type=int, default=8)
    ap.add_argument("--max-pair-bridge-edges", type=int, default=6)
    ap.add_argument("--max-pair-span", type=int, default=320)
    ap.add_argument("--dominance", type=float, default=0.70)
    ap.add_argument("--min-direct", type=int, default=4)
    ap.add_argument("--min-length", type=int, default=200)
    args = ap.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    graph = gp.Graph.from_gfa(args.gfa)
    index = gp.KmerIndex(graph, args.anchor_k)
    base_paths = load_paths(args.base_paths)
    membership = gp.preliminary_membership(base_paths)

    raw_ctx, raw_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, None, args.max_context
    )
    proj_ctx, high_ctx, projection_stats = collect_projection_contexts(
        graph, index, args.projection, args.highk_gfa, out, args.max_context
    )
    second_ctx, second_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, membership, args.max_context
    )
    for key in list(second_ctx):
        baseline = raw_ctx.get(key, 0)
        if second_ctx[key] <= baseline:
            del second_ctx[key]
        else:
            second_ctx[key] -= baseline

    pair_ctx, pair_stats = collect_pair_contexts(
        graph,
        index,
        args.read1,
        args.read2,
        membership,
        args.max_context,
        args.max_pair_bridge_edges,
        args.max_pair_span,
    )
    repeat_ctx = combined_contexts(second_ctx, pair_ctx)
    repeat_paths, repeat_resolve = gp.resolve_paths(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        args.dominance,
        args.min_direct,
        args.min_length,
    )
    repeat_write = gp.write_paths(
        repeat_paths,
        graph,
        out / "stage5_repeat_traversal.fasta",
        out / "stage5_repeat_traversal.paths.tsv",
        args.min_length,
    )

    all_ctx = combined_contexts(raw_ctx, proj_ctx, high_ctx, repeat_ctx)
    simplified, simplify_stats = simplify_graph(graph, all_ctx)
    simplified_paths, simplified_resolve = gp.resolve_paths(
        simplified,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        args.dominance,
        args.min_direct,
        args.min_length,
    )
    simplified_write = gp.write_paths(
        simplified_paths,
        simplified,
        out / "stage6_graph_simplified.fasta",
        out / "stage6_graph_simplified.paths.tsv",
        args.min_length,
    )

    stats = {
        "full_read_threading": raw_stats,
        "projection_threading": projection_stats,
        "second_pass_threading": second_stats,
        "pair_repeat_threading": pair_stats,
        "stage5_repeat_traversal": {**repeat_resolve, **repeat_write},
        "graph_simplification": simplify_stats,
        "stage6_graph_simplified": {**simplified_resolve, **simplified_write},
    }
    (out / "repeat_optimizer_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
