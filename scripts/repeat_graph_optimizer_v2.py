#!/usr/bin/env python3
"""Repeat optimizer v2: context-first traversal from unique flanks.

The baseline graph phaser seeds by coverage*length. Repeats can therefore be
claimed before a path has enough history to phase the correct branch. This
variant keeps the same edge-scoring/safety rules but starts from low-ambiguity
flanks first, allowing full-read/mate/high-k context to be present when a repeat
junction is reached. Nodes are still used at most once; this stage does not
introduce repeat copy-number expansion or new graph edges.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import graph_path_phaser as gp
import repeat_graph_optimizer as rg


def resolve_context_seeded_paths(
    graph: gp.Graph,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    second_ctx: Counter[tuple[str, ...]],
    dominance: float,
    min_direct: int,
    min_length: int,
) -> tuple[list[list[str]], dict[str, int]]:
    used: set[str] = set()
    paths: list[list[str]] = []
    phased_extensions = 0
    branch_stops = 0
    seeds = sorted(
        graph.seqs,
        key=lambda uid: (
            int(graph.ambiguous(uid)),
            max(0, len(graph.inc.get(uid, [])) - 1)
            + max(0, len(graph.out.get(uid, [])) - 1),
            -(graph.coverage.get(uid, 0.0) * len(graph.seqs[uid])),
            -len(graph.seqs[uid]),
            uid,
        ),
    )
    ambiguous_seed_skips = 0
    for seed in seeds:
        if seed in used or graph.rev.get(seed, seed) in used:
            continue
        if graph.ambiguous(seed):
            ambiguous_seed_skips += 1
        path = [seed]
        local_seen = {seed, graph.rev.get(seed, seed)}
        while True:
            current = path[0]
            choice = gp.choose_extension(
                graph,
                path,
                graph.inc.get(current, []),
                used | local_seen,
                False,
                raw_ctx,
                proj_ctx,
                high_ctx,
                second_ctx,
                dominance,
                min_direct,
            )
            if choice is None:
                if len(graph.inc.get(current, [])) > 1:
                    branch_stops += 1
                break
            path.insert(0, choice.uid)
            local_seen.add(choice.uid)
            local_seen.add(graph.rev.get(choice.uid, choice.uid))
            if gp.strong_context(choice):
                phased_extensions += 1
        while True:
            current = path[-1]
            choice = gp.choose_extension(
                graph,
                path,
                graph.out.get(current, []),
                used | local_seen,
                True,
                raw_ctx,
                proj_ctx,
                high_ctx,
                second_ctx,
                dominance,
                min_direct,
            )
            if choice is None:
                if len(graph.out.get(current, [])) > 1:
                    branch_stops += 1
                break
            path.append(choice.uid)
            local_seen.add(choice.uid)
            local_seen.add(graph.rev.get(choice.uid, choice.uid))
            if gp.strong_context(choice):
                phased_extensions += 1
        gp.claim(path, graph, used)
        if len(gp.path_sequence(path, graph)) >= min_length:
            paths.append(path)
    return paths, {
        "paths": len(paths),
        "phased_extensions": phased_extensions,
        "branch_stops": branch_stops,
        "claimed_orientations": len(used),
        "ambiguous_paths_seeded_after_unique": ambiguous_seed_skips,
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
    base_paths = rg.load_paths(args.base_paths)
    membership = gp.preliminary_membership(base_paths)

    raw_ctx, raw_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, None, args.max_context
    )
    proj_ctx, high_ctx, projection_stats = rg.collect_projection_contexts(
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

    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph,
        index,
        args.read1,
        args.read2,
        membership,
        args.max_context,
        args.max_pair_bridge_edges,
        args.max_pair_span,
    )
    repeat_ctx = rg.combined_contexts(second_ctx, pair_ctx)
    repeat_paths, repeat_resolve = resolve_context_seeded_paths(
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
        out / "stage5_repeat_seeded.fasta",
        out / "stage5_repeat_seeded.paths.tsv",
        args.min_length,
    )

    all_ctx = rg.combined_contexts(raw_ctx, proj_ctx, high_ctx, repeat_ctx)
    simplified, simplify_stats = rg.simplify_graph(graph, all_ctx)
    simplified_paths, simplified_resolve = resolve_context_seeded_paths(
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
        out / "stage6_graph_simplified_seeded.fasta",
        out / "stage6_graph_simplified_seeded.paths.tsv",
        args.min_length,
    )

    stats = {
        "full_read_threading": raw_stats,
        "projection_threading": projection_stats,
        "second_pass_threading": second_stats,
        "pair_repeat_threading": pair_stats,
        "stage5_repeat_seeded": {**repeat_resolve, **repeat_write},
        "graph_simplification": simplify_stats,
        "stage6_graph_simplified_seeded": {
            **simplified_resolve,
            **simplified_write,
        },
    }
    (out / "repeat_optimizer_v2_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
