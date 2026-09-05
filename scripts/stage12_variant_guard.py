#!/usr/bin/env python3
"""Stage 12: haplotype-linkage guard for aggressive graph traversal.

Stage11 deliberately explores deeper/looser abundance-flow traversal. Stage12
adds a long-range read-pair linkage layer that can (a) veto a locally attractive
branch contradicted by upstream haplotype markers and (b) rescue a branch that
Stage11 leaves unresolved when one locally eligible child is strongly linked to
the preceding path history.

The linkage signal excludes the current branch source. Therefore it cannot be a
simple re-count of the local PE edge tag: support must connect a candidate child
to one of the previous path unitig families through family-unique exact 21-mers
observed on the same physical read pair.
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from collections import Counter, defaultdict
from pathlib import Path

import graph_path_phaser as gp
import low_abundance_rescue as lr
import repeat_graph_optimizer as rg
import stage11_aggressive_rescue as s11
import stage789_optimizer as s78
import variant_linkage as vl


def choose_extension_variant_guard(
    graph: gp.Graph,
    linkage: vl.MarkerLinkage,
    history: list[str],
    candidates: list[str],
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
) -> tuple[gp.Choice | None, str, list[vl.LinkageScore]]:
    available = [
        uid
        for uid in candidates
        if uid not in used and graph.rev.get(uid, uid) not in used
    ]
    if not available:
        return None, "stop", []

    proposed, stage11_mode = s11.choose_extension_abundance_flow(
        graph,
        history,
        available,
        used,
        forward,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
    )
    decision, _guarded_uid, scores = vl.linkage_decision(
        graph,
        linkage,
        history,
        available,
        forward,
        proposed.uid if proposed is not None else None,
        min_total=3.0,
        min_rescue_support=3.0,
        rescue_share=0.67,
        rescue_margin=1.50,
        veto_share=0.25,
        veto_margin=2.00,
    )
    if decision == "veto":
        return None, "variant_veto", scores
    if proposed is not None:
        confirmed = next((item for item in scores if item.uid == proposed.uid), None)
        if confirmed is not None and confirmed.support >= 3.0 and confirmed.share >= 0.55:
            return proposed, "variant_confirmed", scores
        return proposed, stage11_mode, scores

    local = s11.rank_flow_candidates(
        graph,
        history,
        available,
        used,
        forward,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
    )
    locally_eligible = [item.choice.uid for item in local]
    if not locally_eligible:
        return None, "stop", scores
    decision2, rescued_uid, scores2 = vl.linkage_decision(
        graph,
        linkage,
        history,
        locally_eligible,
        forward,
        None,
        min_total=3.0,
        min_rescue_support=3.0,
        rescue_share=0.67,
        rescue_margin=1.50,
    )
    if decision2 != "rescue" or rescued_uid is None:
        return None, "stop", scores2
    choice = next((item.choice for item in local if item.choice.uid == rescued_uid), None)
    return choice, "variant_rescue", scores2


def resolve_variant_guard_paths(
    graph: gp.Graph,
    linkage: vl.MarkerLinkage,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
    min_length: int,
) -> tuple[list[list[str]], dict[str, int]]:
    used: set[str] = set()
    paths: list[list[str]] = []
    stats: defaultdict[str, int] = defaultdict(int)
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
    for seed in seeds:
        if seed in used or graph.rev.get(seed, seed) in used:
            continue
        path = [seed]
        local_seen = {seed, graph.rev.get(seed, seed)}
        for forward in (False, True):
            while True:
                current = path[-1] if forward else path[0]
                adjacent = graph.out.get(current, []) if forward else graph.inc.get(current, [])
                choice, mode, scores = choose_extension_variant_guard(
                    graph,
                    linkage,
                    path,
                    adjacent,
                    used | local_seen,
                    forward,
                    raw_ctx,
                    proj_ctx,
                    high_ctx,
                    repeat_ctx,
                )
                if mode == "variant_veto":
                    stats["variant_vetoes"] += 1
                if choice is None:
                    if len(adjacent) > 1:
                        stats["branch_stops"] += 1
                    break
                if forward:
                    path.append(choice.uid)
                else:
                    path.insert(0, choice.uid)
                local_seen.update((choice.uid, graph.rev.get(choice.uid, choice.uid)))
                stats[f"{mode}_extensions"] += 1
                if scores and max(item.support for item in scores) >= 3.0:
                    stats["linkage_informative_extensions"] += 1
        gp.claim(path, graph, used)
        if len(gp.path_sequence(path, graph)) >= min_length:
            paths.append(path)
    stats["paths"] = len(paths)
    stats["claimed_orientations"] = len(used)
    return paths, dict(stats)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--segment-anchor-bases", type=int, default=31)
    ap.add_argument("--marker-k", type=int, default=21)
    args = ap.parse_args()

    started = time.monotonic()
    scripts = Path(__file__).resolve().parent
    out = args.pipeline_dir
    base = out / "current_pipeline"
    graph_opt = out / "graph_optimizer"
    repeat_opt = out / "repeat_optimizer"
    current = base / "step6_strain_projection.fasta"
    backbone = out / "bridge_backbone.fasta"
    target_gfa = base / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = base / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"
    base_paths = graph_opt / "stage4_second_pass.paths.tsv"
    consensus_add = out / "stage11_aggressive" / "rare_consensus_additions.fasta"
    required = [
        current,
        backbone,
        target_gfa,
        projection_primary,
        highk_gfa,
        base_paths,
        consensus_add,
        args.read1,
        args.read2,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage12 inputs: " + ", ".join(missing))

    stage12 = out / "stage12_variant_guard"
    stage12.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    graph = gp.Graph.from_gfa(target_gfa)
    index = gp.KmerIndex(graph, 31)
    membership = gp.preliminary_membership(rg.load_paths(base_paths))
    raw_ctx, raw_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, None, 10
    )
    proj_ctx, high_ctx, projection_stats = rg.collect_projection_contexts(
        graph,
        index,
        [projection_primary, projection_haplotigs],
        [highk_gfa],
        repeat_opt,
        10,
    )
    second_ctx, second_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, membership, 10
    )
    for key in list(second_ctx):
        baseline = raw_ctx.get(key, 0)
        if second_ctx[key] <= baseline:
            del second_ctx[key]
        else:
            second_ctx[key] -= baseline
    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph, index, args.read1, args.read2, membership, 10, 8, 420
    )
    repeat_ctx = rg.combined_contexts(second_ctx, pair_ctx)
    all_ctx = rg.combined_contexts(raw_ctx, proj_ctx, high_ctx, repeat_ctx)
    simplified, simplify_stats = rg.simplify_graph(graph, all_ctx)

    linkage_started = time.monotonic()
    linkage = vl.collect_linkage(
        simplified,
        args.read1,
        args.read2,
        k=args.marker_k,
        stride=2,
        min_markers=2,
        max_families=8,
    )
    timings["variant_linkage_collection"] = time.monotonic() - linkage_started

    paths, path_stats = resolve_variant_guard_paths(
        simplified,
        linkage,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        200,
    )
    raw = stage12 / "variant_guard_flow.raw.fasta"
    write_stats = gp.write_paths(
        paths,
        simplified,
        raw,
        stage12 / "variant_guard_flow.paths.tsv",
        200,
    )
    flow = s78.emit_stage(
        scripts,
        raw,
        current,
        out,
        "stage12_variant_guard_flow",
        args.segment_anchor_bases,
        timings,
    )
    flow_rare = lr.make_union_candidate(
        scripts,
        flow,
        [consensus_add],
        stage12 / "candidate_variant_flow_rare",
        timings,
    )
    flow_local = s11.positive_gap_localfill(
        scripts,
        flow,
        args.read1,
        args.read2,
        stage12 / "candidate_variant_flow_localfill",
        args.threads,
        timings,
    )
    flow_rare_local = s11.positive_gap_localfill(
        scripts,
        flow_rare,
        args.read1,
        args.read2,
        stage12 / "candidate_variant_flow_rare_localfill",
        args.threads,
        timings,
    )

    stats = {
        "pipeline": "bridge-stage12-variant-guard-v1",
        "policy": {
            "production_backbone_replaced": False,
            "local_edge_required_for_variant_rescue": True,
            "branch_source_excluded_from_linkage": True,
            "marker_k": args.marker_k,
            "min_unique_markers_per_family_fragment": 2,
        },
        "linkage": {
            "fragments": linkage.fragments,
            "informative_fragments": linkage.informative_fragments,
            "family_pairs": len(linkage.pair_counts),
            "families": len(linkage.family_depth),
        },
        "path_stats": path_stats,
        "write_stats": write_stats,
        "simplification": simplify_stats,
        "contexts": {
            "raw": raw_stats,
            "projection": projection_stats,
            "second": second_stats,
            "pair": pair_stats,
        },
        "outputs": {
            "variant_flow": str(flow),
            "variant_flow_localfill": str(flow_local),
            "variant_flow_rare": str(flow_rare),
            "variant_flow_rare_localfill": str(flow_rare_local),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (stage12 / "stage12_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
