#!/usr/bin/env python3
"""Stage20: resolve only the gaps between trustworthy exact k31 anchors.

Stage16 showed that free k15 chaining is too ambiguous. Stage17 preserves exact
threading, but its fallback is used only when a read has no exact multi-unitig
segment. A common remaining failure mode is a read that has exact k31 anchors on
both sides of a small ambiguous graph bubble: exact threading splits at the
bubble, while unconstrained low-k chaining searches far too much of the graph.

This module keeps the exact anchors fixed and enumerates only bounded source to
sink graph paths between consecutive exact anchors. Candidate paths are ranked
by discriminative k19 evidence from the read interval; k15 is used only when
k19 has no decisive evidence. Flank k-mers are removed before scoring, so low-k
evidence must come from the alternative internal path itself. No graph edge is
invented. Accepted contexts are add-only on top of exact read threading and
final sequence is add-only on top of the Stage10 strict baseline.

No reference is used.
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import graph_path_phaser as gp
import repeat_graph_optimizer as rg
import stage14_amplified_methods as s14
import stage16_root_cause as s16
import stage17_frontier_recovery as s17
import stage789_optimizer as s78


@dataclass(frozen=True)
class ExactAnchor:
    start_pos: int
    end_pos: int
    uid: str


@dataclass(frozen=True)
class GapResolution:
    source: str
    target: str
    source_pos: int
    target_pos: int
    path: tuple[str, ...]
    evidence_k: int
    best_hits: int
    second_hits: int
    physical_edges: int
    total_edges: int


def exact_anchor_runs(seq: str, index: gp.KmerIndex) -> list[ExactAnchor]:
    """Return start/end positions for each consecutive unique-unitig run."""
    anchors: list[ExactAnchor] = []
    for pos, key in gp.rolling_keys(seq, index.k):
        uid = index.unique.get(key)
        if uid is None:
            continue
        if anchors and anchors[-1].uid == uid:
            old = anchors[-1]
            anchors[-1] = ExactAnchor(old.start_pos, pos, uid)
        else:
            anchors.append(ExactAnchor(pos, pos, uid))
    return anchors


def enumerate_bounded_paths(
    graph: gp.Graph,
    source: str,
    target: str,
    *,
    max_edges: int = 5,
    max_paths: int = 8,
    max_states: int = 160,
    max_internal_bp: int = 220,
) -> list[list[str]]:
    """Enumerate a small acyclic source->target path family."""
    if source == target:
        return [[source]]
    queue: list[tuple[str, list[str], int]] = [(source, [source], 0)]
    out: list[list[str]] = []
    states = 0
    while queue and states < max_states and len(out) < max_paths:
        node, path, internal_bp = queue.pop(0)
        states += 1
        edges = len(path) - 1
        if edges >= max_edges:
            continue
        for child in graph.out.get(node, []):
            if child in path or graph.rev.get(child, child) in path:
                continue
            path2 = path + [child]
            if child == target:
                out.append(path2)
                if len(out) >= max_paths:
                    break
                continue
            bp2 = internal_bp + max(1, len(graph.seqs[child]) - graph.k)
            if bp2 > max_internal_bp:
                continue
            queue.append((child, path2, bp2))
    return out


def keyset(seq: str, k: int) -> set[int]:
    return {key for _pos, key in gp.rolling_keys(seq, k)}


def internal_path_keys(path: list[str], graph: gp.Graph, k: int) -> set[int]:
    keys: set[int] = set()
    for uid in path[1:-1]:
        keys.update(keyset(graph.seqs[uid], k))
    return keys


def physical_edge_count(path: list[str], graph: gp.Graph) -> tuple[int, int]:
    supported = 0
    total = 0
    for left, right in zip(path, path[1:]):
        total += 1
        ev = graph.edge.get((left, right), gp.EdgeEvidence())
        supported += int((ev.direct + ev.gapped + ev.pairs) > 0)
    return supported, total


def choose_path_by_read(
    seq: str,
    left: ExactAnchor,
    right: ExactAnchor,
    paths: list[list[str]],
    graph: gp.Graph,
) -> GapResolution | None:
    if len(paths) < 2:
        return None
    # Constrain the scoring interval to the actual unresolved gap: the terminal
    # exact k31 of the left run to the first exact k31 of the right run.
    start = max(0, left.end_pos - 8)
    end = min(len(seq), right.start_pos + graph.k + 8)
    if end <= start:
        return None
    segment = seq[start:end]

    for k, min_hits, margin in ((19, 3, 2), (15, 6, 4)):
        read_keys = keyset(segment, k)
        flank_keys = keyset(graph.seqs[left.uid], k) | keyset(
            graph.seqs[right.uid], k
        )
        discriminative = read_keys - flank_keys
        if not discriminative:
            continue
        scores: list[tuple[int, int, int, int, tuple[str, ...], list[str]]] = []
        for path in paths:
            hits = len(discriminative & internal_path_keys(path, graph, k))
            physical, total = physical_edge_count(path, graph)
            scores.append((hits, physical, -len(path), total, tuple(path), path))
        scores.sort(reverse=True)
        best = scores[0]
        second_hits = scores[1][0]
        best_hits, physical, _neg_len, total, _sig, best_path = best
        decisive = best_hits >= min_hits and best_hits >= second_hits + margin
        if not decisive:
            continue
        if total and physical / total < 0.75:
            continue
        return GapResolution(
            left.uid,
            right.uid,
            left.end_pos,
            right.start_pos,
            tuple(best_path),
            k,
            best_hits,
            second_hits,
            physical,
            total,
        )
    return None


def resolve_read_anchor_gaps(
    seq: str,
    graph: gp.Graph,
    exact_index: gp.KmerIndex,
    *,
    max_edges: int = 5,
) -> list[GapResolution]:
    anchors = exact_anchor_runs(seq, exact_index)
    resolved: list[GapResolution] = []
    for left, right in zip(anchors, anchors[1:]):
        if left.uid == right.uid or right.start_pos <= left.end_pos:
            continue
        if right.uid in graph.out.get(left.uid, []):
            continue
        if gp.unique_short_bridge(graph, left.uid, right.uid, 3) is not None:
            continue
        read_delta = right.start_pos - left.end_pos
        paths = enumerate_bounded_paths(
            graph,
            left.uid,
            right.uid,
            max_edges=max_edges,
            max_internal_bp=read_delta + 2 * graph.k + 20,
        )
        if len(paths) < 2:
            continue
        choice = choose_path_by_read(seq, left, right, paths, graph)
        if choice is not None:
            resolved.append(choice)
    return resolved


def collect_anchor_gap_contexts(
    graph: gp.Graph,
    exact_index: gp.KmerIndex,
    read1: Path,
    read2: Path,
    *,
    max_context: int = 10,
) -> tuple[Counter[tuple[str, ...]], dict[str, object], list[GapResolution]]:
    contexts: Counter[tuple[str, ...]] = Counter()
    all_resolutions: list[GapResolution] = []
    stats: dict[str, object] = {
        "reads": 0,
        "reads_with_exact_anchor_runs": 0,
        "reads_with_resolved_gap": 0,
        "resolved_gaps": 0,
        "resolved_k19": 0,
        "resolved_k15": 0,
    }
    path_counts: Counter[tuple[str, ...]] = Counter()
    for fastq in (read1, read2):
        for _name, seq in gp.read_fastq(fastq):
            stats["reads"] = int(stats["reads"]) + 1
            anchors = exact_anchor_runs(seq, exact_index)
            if len(anchors) >= 2:
                stats["reads_with_exact_anchor_runs"] = int(
                    stats["reads_with_exact_anchor_runs"]
                ) + 1
            resolutions = resolve_read_anchor_gaps(seq, graph, exact_index)
            if not resolutions:
                continue
            stats["reads_with_resolved_gap"] = int(stats["reads_with_resolved_gap"]) + 1
            for item in resolutions:
                stats["resolved_gaps"] = int(stats["resolved_gaps"]) + 1
                stats[f"resolved_k{item.evidence_k}"] = int(
                    stats[f"resolved_k{item.evidence_k}"]
                ) + 1
                path_counts[item.path] += 1
                all_resolutions.append(item)
                gp.add_context(contexts, list(item.path), max_context)
    stats["distinct_resolved_paths"] = len(path_counts)
    stats["paths_supported_by_2plus_reads"] = sum(v >= 2 for v in path_counts.values())
    stats["contexts"] = len(contexts)
    stats["context_weight"] = sum(contexts.values())
    return contexts, stats, all_resolutions


def write_resolutions(items: list[GapResolution], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter(item.path for item in items)
    with path.open("w") as handle:
        handle.write(
            "path_support\tevidence_k\tbest_hits\tsecond_hits\tphysical_edges\t"
            "total_edges\tsource_pos\ttarget_pos\tpath\n"
        )
        for item in sorted(
            items,
            key=lambda x: (
                -counts[x.path],
                -x.evidence_k,
                -x.best_hits,
                x.path,
                x.source_pos,
            ),
        ):
            handle.write(
                f"{counts[item.path]}\t{item.evidence_k}\t{item.best_hits}\t"
                f"{item.second_hits}\t{item.physical_edges}\t{item.total_edges}\t"
                f"{item.source_pos}\t{item.target_pos}\t{','.join(item.path)}\n"
            )


def build_anchor_gap_candidate(
    scripts: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    read1: Path,
    read2: Path,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    base = pipeline_dir / "current_pipeline"
    graph_opt = pipeline_dir / "graph_optimizer"
    repeat_opt = pipeline_dir / "repeat_optimizer"
    target_gfa = base / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = base / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"
    base_paths = graph_opt / "stage4_second_pass.paths.tsv"

    graph = gp.Graph.from_gfa(target_gfa)
    membership = gp.preliminary_membership(rg.load_paths(base_paths))
    exact_index = gp.KmerIndex(graph, 31)

    started = time.monotonic()
    exact_ctx, exact_stats = gp.collect_read_contexts(
        graph, exact_index, read1, read2, None, 10
    )
    gap_ctx, gap_stats, resolutions = collect_anchor_gap_contexts(
        graph, exact_index, read1, read2, max_context=10
    )
    increment = s16.conservative_increment(exact_ctx, gap_ctx, max_weight=6)
    enhanced_raw = Counter(exact_ctx)
    enhanced_raw.update(increment)
    timings["anchor_gap_threading"] = time.monotonic() - started

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
    all_ctx = rg.combined_contexts(enhanced_raw, proj_ctx, high_ctx, repeat_ctx)
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

    outdir = pipeline_dir / "stage20_anchor_gap_thread"
    outdir.mkdir(parents=True, exist_ok=True)
    write_resolutions(resolutions, outdir / "resolved_anchor_gaps.tsv")
    selected = s17.select_hybrid_path_additions(
        paths, simplified, strict_baseline, gap_ctx
    )
    s17.write_hybrid_evidence(selected, outdir / "anchor_gap_selected.tsv")
    additions = outdir / "anchor_gap_additions.fasta"
    s14.write_fasta(
        ((f"stage20_anchor_gap_{idx:06d}", item.seq) for idx, item in enumerate(selected, 1)),
        additions,
    )
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_anchor_gap",
        timings,
        min_overlap=81,
    )
    return final, additions, {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "exact_threading": exact_stats,
        "anchor_gap": gap_stats,
        "incremental_contexts": len(increment),
        "incremental_context_weight": sum(increment.values()),
        "projection": projection_stats,
        "second_pass": second_stats,
        "pair_threading": pair_stats,
        "simplification": simplify_stats,
        "path_resolution": resolve_stats,
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "selected_fresh31": sum(item.fresh31 for item in selected),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    args = ap.parse_args()

    started = time.monotonic()
    scripts = Path(__file__).resolve().parent
    stage10 = args.pipeline_dir / "stage10_multik_rescue"
    strict_baseline = stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    required = [
        args.read1,
        args.read2,
        strict_baseline,
        args.pipeline_dir / "current_pipeline" / "iterative" / "k31_resolve" / "assembly.gfa",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage20 inputs: " + ", ".join(missing))

    timings: dict[str, float] = {}
    final, additions, method = build_anchor_gap_candidate(
        scripts,
        args.pipeline_dir,
        strict_baseline,
        args.read1,
        args.read2,
        timings,
    )
    root = args.pipeline_dir / "stage20_anchor_gap_thread"
    stats = {
        "pipeline": "bridge-stage20-anchor-gap-thread-v2",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "metric_targets": False,
            "fixed_evidence": "unique exact k31 anchor runs",
            "search_space": "bounded graph paths from left-run terminal anchor to right-run initial anchor",
            "branch_evidence": "flank-subtracted k19, k15 fallback",
            "graph_edges": "existing GFA edges only",
            "output": "add-only on Stage10 strict",
        },
        "method": method,
        "outputs": {"final": str(final), "additions": str(additions)},
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (root / "stage20_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
