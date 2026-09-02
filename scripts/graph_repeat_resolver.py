#!/usr/bin/env python3
"""Repeat-aware graph traversal after full-read phasing.

Two cumulative stages are emitted:
  5. bounded lookahead traversal across ambiguous repeat/bubble junctions
  6. conservative simplification of short reconvergent bubbles, followed by
     the same lookahead traversal

No new graph edge is invented. Lookahead only follows existing GFA edges and
requires context spanning beyond the immediate branch. Simplification only
blocks losing first edges of short reconvergent bubbles with a strong evidence
margin, so low-coverage disconnected sequence is not globally pruned.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from graph_path_phaser import (
    EdgeEvidence,
    Graph,
    KmerIndex,
    collect_fasta_contexts,
    collect_read_contexts,
    preliminary_membership,
    write_paths,
)


def read_paths(path: Path) -> list[list[str]]:
    paths: list[list[str]] = []
    with path.open() as handle:
        header = next(handle, None)
        if header is None:
            return paths
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            nodes = [x for x in fields[3].split(",") if x]
            if nodes:
                paths.append(nodes)
    return paths


def crossing_context(
    counter: Counter[tuple[str, ...]],
    history: list[str],
    extension: list[str],
    forward: bool,
) -> tuple[int, int]:
    best_support = 0
    best_len = 0
    if not history or not extension:
        return 0, 0
    if forward:
        left = history[-5:]
        right = extension[:5]
    else:
        left = extension[-5:]
        right = history[:5]
    for l_take in range(1, min(5, len(left)) + 1):
        for r_take in range(1, min(5, len(right)) + 1):
            key = tuple(left[-l_take:] + right[:r_take])
            if len(key) < 3:
                continue
            support = counter.get(key, 0)
            if support > 0 and (
                len(key) > best_len
                or (len(key) == best_len and support > best_support)
            ):
                best_support = support
                best_len = len(key)
    return best_support, best_len


def enumerate_extensions(
    graph: Graph,
    start: str,
    forward: bool,
    max_edges: int,
    forbidden: set[str],
    blocked: set[tuple[str, str]],
) -> list[list[str]]:
    neighbors = graph.out if forward else graph.inc
    frontier: list[list[str]] = [
        [uid]
        for uid in neighbors[start]
        if uid not in forbidden
        and ((start, uid) if forward else (uid, start)) not in blocked
    ]
    out: list[list[str]] = []
    for _ in range(max_edges):
        if not frontier:
            break
        nxt: list[list[str]] = []
        for ext in frontier:
            out.append(ext)
            node = ext[-1]
            for child in neighbors[node]:
                edge = (node, child) if forward else (child, node)
                if edge in blocked or child in forbidden or child in ext:
                    continue
                nxt.append(ext + [child])
        frontier = nxt
    return out


def extension_score(
    graph: Graph,
    history: list[str],
    ext: list[str],
    forward: bool,
    raw_ctx,
    proj_ctx,
    high_ctx,
    second_ctx,
) -> tuple[float, tuple[int, int, int, int, int, int, int, int]]:
    # enumerate_extensions stores reverse traversal as near->far; convert it to
    # genomic path order far->near before scoring contexts and graph edges.
    oriented_ext = ext if forward else list(reversed(ext))
    rs, rl = crossing_context(raw_ctx, history, oriented_ext, forward)
    ps, pl = crossing_context(proj_ctx, history, oriented_ext, forward)
    hs, hl = crossing_context(high_ctx, history, oriented_ext, forward)
    ss, sl = crossing_context(second_ctx, history, oriented_ext, forward)
    edge_score = 0.0
    chain = (
        [history[-1]] + oriented_ext
        if forward
        else oriented_ext + [history[0]]
    )
    for a, b in zip(chain, chain[1:]):
        ev = graph.edge.get((a, b), EdgeEvidence())
        edge_score += ev.direct * 1.5 + ev.gapped * 0.5 + ev.pairs * 2.0
    score = (
        edge_score
        + rs * max(3, rl) * 9.0
        + ps * max(3, pl) * 6.0
        + hs * max(3, hl) * 11.0
        + ss * max(3, sl) * 12.0
    )
    return score, (rs, rl, ps, pl, hs, hl, ss, sl)


def choose_lookahead(
    graph: Graph,
    history: list[str],
    forward: bool,
    forbidden: set[str],
    blocked: set[tuple[str, str]],
    raw_ctx,
    proj_ctx,
    high_ctx,
    second_ctx,
    max_edges: int,
    margin: float,
) -> list[str] | None:
    start = history[-1] if forward else history[0]
    immediate = graph.out[start] if forward else graph.inc[start]
    available = [
        u
        for u in immediate
        if u not in forbidden
        and ((start, u) if forward else (u, start)) not in blocked
    ]
    if len(available) < 2:
        return None
    scored = []
    for ext in enumerate_extensions(
        graph, start, forward, max_edges, forbidden, blocked
    ):
        score, ctx = extension_score(
            graph,
            history,
            ext,
            forward,
            raw_ctx,
            proj_ctx,
            high_ctx,
            second_ctx,
        )
        rs, rl, ps, pl, hs, hl, ss, sl = ctx
        spanning = (
            (rl >= 3 and rs >= 2)
            or (sl >= 3 and ss >= 2)
            or (hl >= 3 and hs >= 1)
            or (pl >= 3 and ps >= 1)
        )
        if spanning:
            scored.append((score, ext, ctx))
    if not scored:
        return None
    scored.sort(key=lambda item: (-item[0], -len(item[1]), item[1]))
    best_score, best_ext, _ = scored[0]
    competing = [item for item in scored[1:] if item[1][0] != best_ext[0]]
    second = max((item[0] for item in competing), default=0.0)
    if best_score <= 0 or (second > 0 and best_score < second * margin):
        return None
    same = [
        item
        for item in scored
        if item[1][0] == best_ext[0] and item[0] >= best_score * 0.70
    ]
    same.sort(key=lambda item: (-len(item[1]), -item[0], item[1]))
    return same[0][1]


def edge_choice(
    graph: Graph,
    history: list[str],
    forward: bool,
    forbidden: set[str],
    blocked: set[tuple[str, str]],
    raw_ctx,
    proj_ctx,
    high_ctx,
    second_ctx,
    dominance: float,
    min_direct: int,
):
    start = history[-1] if forward else history[0]
    neighbors = graph.out[start] if forward else graph.inc[start]
    choices = []
    for uid in neighbors:
        edge = (start, uid) if forward else (uid, start)
        if uid in forbidden or edge in blocked:
            continue
        score, ctx = extension_score(
            graph,
            history,
            [uid],
            forward,
            raw_ctx,
            proj_ctx,
            high_ctx,
            second_ctx,
        )
        ev = graph.edge.get(edge, EdgeEvidence())
        choices.append((score, uid, ctx, ev.direct))
    if not choices:
        return None
    if len(choices) == 1:
        return [choices[0][1]] if choices[0][0] > 0 or len(neighbors) == 1 else None
    choices.sort(key=lambda x: (-x[0], x[1]))
    best = choices[0]
    total = sum(max(0.0, c[0]) for c in choices)
    frac = best[0] / total if total > 0 else 0.0
    rs, rl, ps, pl, hs, hl, ss, sl = best[2]
    strong = (
        (rl >= 3 and rs >= 2)
        or (sl >= 3 and ss >= 2)
        or (hl >= 3 and hs >= 1)
        or (pl >= 3 and ps >= 1)
    )
    if strong and (len(choices) == 1 or best[0] > choices[1][0] * 1.10):
        return [best[1]]
    if best[3] >= min_direct and frac >= dominance:
        return [best[1]]
    return None


def resolve(
    graph: Graph,
    raw_ctx,
    proj_ctx,
    high_ctx,
    second_ctx,
    blocked: set[tuple[str, str]],
    dominance: float,
    min_direct: int,
    min_length: int,
    lookahead_edges: int,
    lookahead_margin: float,
):
    used: set[str] = set()
    paths: list[list[str]] = []
    branch_stops = lookahead_extensions = 0
    seeds = sorted(
        graph.seqs,
        key=lambda uid: (
            -(graph.coverage.get(uid, 0.0) * len(graph.seqs[uid])),
            -len(graph.seqs[uid]),
            uid,
        ),
    )
    for seed in seeds:
        if seed in used or graph.rev.get(seed, seed) in used:
            continue
        path = [seed]
        local = {seed, graph.rev.get(seed, seed)}
        for forward in (False, True):
            while True:
                start = path[-1] if forward else path[0]
                neighbors = graph.out[start] if forward else graph.inc[start]
                forbidden = used | local
                ext = edge_choice(
                    graph,
                    path,
                    forward,
                    forbidden,
                    blocked,
                    raw_ctx,
                    proj_ctx,
                    high_ctx,
                    second_ctx,
                    dominance,
                    min_direct,
                )
                if ext is None and len(neighbors) > 1:
                    ext = choose_lookahead(
                        graph,
                        path,
                        forward,
                        forbidden,
                        blocked,
                        raw_ctx,
                        proj_ctx,
                        high_ctx,
                        second_ctx,
                        lookahead_edges,
                        lookahead_margin,
                    )
                    if ext is not None:
                        lookahead_extensions += 1
                if ext is None:
                    if len(neighbors) > 1:
                        branch_stops += 1
                    break
                ordered = ext if forward else list(reversed(ext))
                if forward:
                    path.extend(ordered)
                else:
                    path[:0] = ordered
                for uid in ordered:
                    local.add(uid)
                    local.add(graph.rev.get(uid, uid))
        for uid in path:
            used.add(uid)
            used.add(graph.rev.get(uid, uid))
        seq_len = len(graph.seqs[path[0]]) + sum(
            max(0, len(graph.seqs[u]) - graph.k) for u in path[1:]
        )
        if seq_len >= min_length:
            paths.append(path)
    return paths, {
        "paths": len(paths),
        "branch_stops": branch_stops,
        "lookahead_extensions": lookahead_extensions,
        "blocked_edges": len(blocked),
        "claimed_orientations": len(used),
    }


def paths_from_branch(
    graph: Graph, source: str, max_edges: int, max_bp: int
) -> dict[str, list[list[str]]]:
    by_sink: dict[str, list[list[str]]] = defaultdict(list)
    frontier = [[source, child] for child in graph.out[source]]
    for _ in range(max_edges):
        nxt = []
        for path in frontier:
            node = path[-1]
            bp = len(graph.seqs[path[0]]) + sum(
                max(0, len(graph.seqs[u]) - graph.k) for u in path[1:]
            )
            if bp > max_bp:
                continue
            if node != source and len(graph.inc[node]) > 1:
                by_sink[node].append(path)
            for child in graph.out[node]:
                if child not in path:
                    nxt.append(path + [child])
        frontier = nxt
    return by_sink


def simplify_bubbles(
    graph: Graph,
    raw_ctx,
    proj_ctx,
    high_ctx,
    second_ctx,
    max_edges: int,
    max_bp: int,
    margin: float,
) -> tuple[set[tuple[str, str]], dict[str, int]]:
    blocked: set[tuple[str, str]] = set()
    bubbles = simplified = 0
    for source in graph.seqs:
        if len(graph.out[source]) < 2:
            continue
        by_sink = paths_from_branch(graph, source, max_edges, max_bp)
        for _sink, paths in by_sink.items():
            firsts = {p[1] for p in paths if len(p) >= 2}
            if len(firsts) < 2:
                continue
            bubbles += 1
            candidates = []
            for path in paths:
                ext = path[1:]
                score, ctx = extension_score(
                    graph,
                    [source],
                    ext,
                    True,
                    raw_ctx,
                    proj_ctx,
                    high_ctx,
                    second_ctx,
                )
                first_cov = graph.coverage.get(ext[0], 0.0)
                candidates.append((score, first_cov, ext, ctx))
            candidates.sort(key=lambda x: (-x[0], -x[1], x[2]))
            best = candidates[0]
            alternatives = [x for x in candidates[1:] if x[2][0] != best[2][0]]
            if not alternatives:
                continue
            second = max(alternatives, key=lambda x: x[0])
            rs, rl, ps, pl, hs, hl, ss, sl = best[3]
            strong = (
                (rl >= 3 and rs >= 2)
                or (sl >= 3 and ss >= 2)
                or (hl >= 3 and hs >= 1)
                or (pl >= 3 and ps >= 1)
            )
            if not strong or best[0] <= 0 or best[0] < second[0] * margin:
                continue
            if (
                second[1] > 0
                and best[1] > 0
                and second[1] > best[1] * 0.60
            ):
                continue
            for alt in alternatives:
                if alt[2][0] != best[2][0]:
                    blocked.add((source, alt[2][0]))
            simplified += 1
            break
    return blocked, {
        "candidate_bubbles": bubbles,
        "simplified_bubbles": simplified,
        "blocked_edges": len(blocked),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path)
    ap.add_argument("--projection", type=Path, action="append", default=[])
    ap.add_argument("--highk-gfa", type=Path, action="append", default=[])
    ap.add_argument("--preliminary-paths", type=Path, required=True)
    ap.add_argument("-o", "--output-dir", type=Path, required=True)
    ap.add_argument("--anchor-k", type=int, default=31)
    ap.add_argument("--max-context", type=int, default=6)
    ap.add_argument("--dominance", type=float, default=0.72)
    ap.add_argument("--min-direct", type=int, default=4)
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--lookahead-edges", type=int, default=5)
    ap.add_argument("--lookahead-margin", type=float, default=1.30)
    ap.add_argument("--bubble-edges", type=int, default=4)
    ap.add_argument("--bubble-bp", type=int, default=800)
    ap.add_argument("--bubble-margin", type=float, default=3.0)
    args = ap.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    graph = Graph.from_gfa(args.gfa)
    index = KmerIndex(graph, args.anchor_k)
    raw_ctx, raw_stats = collect_read_contexts(
        graph, index, args.read1, args.read2, None, args.max_context
    )
    proj_ctx, proj_stats = collect_fasta_contexts(
        graph, index, args.projection, args.max_context
    )

    high_ctx: Counter[tuple[str, ...]] = Counter()
    high_stats = {
        "records": 0,
        "segments": 0,
        "accepted_segments": 0,
        "contexts": 0,
    }
    for high_gfa in args.highk_gfa:
        hg = Graph.from_gfa(high_gfa)
        tmp = out / (high_gfa.parent.name + ".highk_unitigs.fasta")
        with tmp.open("w") as handle:
            for uid, seq in hg.seqs.items():
                handle.write(f">{uid}\n{seq}\n")
        ctx, st = collect_fasta_contexts(
            graph, index, [tmp], args.max_context, only_ambiguous=True
        )
        high_ctx.update(ctx)
        for key in ("records", "segments", "accepted_segments"):
            high_stats[key] += st.get(key, 0)
    high_stats["contexts"] = len(high_ctx)

    prelim = read_paths(args.preliminary_paths)
    membership = preliminary_membership(prelim)
    second_ctx, second_stats = collect_read_contexts(
        graph,
        index,
        args.read1,
        args.read2,
        membership,
        args.max_context,
    )
    for key in list(second_ctx):
        base = raw_ctx.get(key, 0)
        if second_ctx[key] <= base:
            del second_ctx[key]
        else:
            second_ctx[key] -= base
    second_stats["novel_contexts"] = len(second_ctx)
    second_stats["novel_support"] = sum(second_ctx.values())

    stage5, st5 = resolve(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        second_ctx,
        set(),
        args.dominance,
        args.min_direct,
        args.min_length,
        args.lookahead_edges,
        args.lookahead_margin,
    )
    st5.update(
        write_paths(
            stage5,
            graph,
            out / "stage5_repeat_traversal.fasta",
            out / "stage5_repeat_traversal.paths.tsv",
            args.min_length,
        )
    )

    blocked, simpl_stats = simplify_bubbles(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        second_ctx,
        args.bubble_edges,
        args.bubble_bp,
        args.bubble_margin,
    )
    stage6, st6 = resolve(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        second_ctx,
        blocked,
        args.dominance,
        args.min_direct,
        args.min_length,
        args.lookahead_edges,
        args.lookahead_margin,
    )
    st6.update(
        write_paths(
            stage6,
            graph,
            out / "stage6_graph_simplified.fasta",
            out / "stage6_graph_simplified.paths.tsv",
            args.min_length,
        )
    )

    stats = {
        "full_read_threading": raw_stats,
        "projection_threading": proj_stats,
        "highk_threading": high_stats,
        "second_pass_threading": second_stats,
        "stage5_repeat_traversal": st5,
        "simplification": simpl_stats,
        "stage6_graph_simplified": st6,
    }
    (out / "repeat_resolver_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
