#!/usr/bin/env python3
"""Promote the proven stage-8 traversal as Bridge's production backbone.

This module intentionally stops after stage 8.  The legacy stage-9 residual
patch experiment remains in stage789_optimizer.py for reproducibility, but is
not part of the promoted pipeline because the 200k Zymo benchmark increased
misassemblies without a meaningful genome-fraction gain.
"""
from __future__ import annotations

import argparse
import json
import resource
import shutil
import time
from collections import Counter
from pathlib import Path

import graph_path_phaser as gp
import repeat_graph_optimizer as rg
import stage789_optimizer as s78


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--segment-anchor-bases", type=int, default=31)
    ap.add_argument("--lookahead-depth", type=int, default=3)
    ap.add_argument("--lookahead-max-branch", type=int, default=4)
    ap.add_argument("--lookahead-discount", type=float, default=0.70)
    ap.add_argument("--lookahead-dominance", type=float, default=0.60)
    ap.add_argument("--lookahead-margin", type=float, default=1.15)
    ap.add_argument("--dominance", type=float, default=0.70)
    ap.add_argument("--min-direct", type=int, default=4)
    ap.add_argument("--min-length", type=int, default=200)
    args = ap.parse_args()

    if args.segment_anchor_bases < 31:
        raise SystemExit("segment-anchor-bases must be >=31")
    if args.lookahead_depth < 2 or args.lookahead_max_branch < 2:
        raise SystemExit("lookahead depth/max-branch must both be >=2")
    if not 0.0 < args.lookahead_discount <= 1.0:
        raise SystemExit("lookahead-discount must be in (0,1]")
    if not 0.5 <= args.lookahead_dominance <= 1.0:
        raise SystemExit("lookahead-dominance must be in [0.5,1]")
    if args.lookahead_margin < 1.0:
        raise SystemExit("lookahead-margin must be >=1")

    scripts = Path(__file__).resolve().parent
    out = args.pipeline_dir
    base = out / "current_pipeline"
    graph_opt = out / "graph_optimizer"
    repeat_opt = out / "repeat_optimizer"
    current = base / "step6_strain_projection.fasta"
    target_gfa = base / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = base / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"
    base_paths = graph_opt / "stage4_second_pass.paths.tsv"
    required = [current, target_gfa, projection_primary, highk_gfa, base_paths]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing stage1-6 pipeline inputs: " + ", ".join(missing))

    started = time.monotonic()
    timings: dict[str, float] = {}
    graph = gp.Graph.from_gfa(target_gfa)
    index = gp.KmerIndex(graph, 31)
    membership = gp.preliminary_membership(rg.load_paths(base_paths))

    raw_ctx, raw_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, None, 8
    )
    proj_ctx, high_ctx, projection_stats = rg.collect_projection_contexts(
        graph,
        index,
        [projection_primary, projection_haplotigs],
        [highk_gfa],
        repeat_opt,
        8,
    )
    second_ctx, second_stats = gp.collect_read_contexts(
        graph, index, args.read1, args.read2, membership, 8
    )
    for key in list(second_ctx):
        baseline = raw_ctx.get(key, 0)
        if second_ctx[key] <= baseline:
            del second_ctx[key]
        else:
            second_ctx[key] -= baseline
    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph, index, args.read1, args.read2, membership, 8, 6, 320
    )
    repeat_ctx = rg.combined_contexts(second_ctx, pair_ctx)
    all_ctx: Counter[tuple[str, ...]] = rg.combined_contexts(
        raw_ctx, proj_ctx, high_ctx, repeat_ctx
    )
    simplified, simplify_stats = rg.simplify_graph(graph, all_ctx)

    stage7_paths, stage7_stats = s78.resolve_lookahead_seeded_paths(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        args.dominance,
        args.min_direct,
        args.min_length,
        args.lookahead_depth,
        args.lookahead_max_branch,
        args.lookahead_discount,
        args.lookahead_dominance,
        args.lookahead_margin,
    )
    stage7_raw = repeat_opt / "stage7_bounded_lookahead.fasta"
    stage7_write = gp.write_paths(
        stage7_paths,
        graph,
        stage7_raw,
        repeat_opt / "stage7_bounded_lookahead.paths.tsv",
        args.min_length,
    )
    stage7 = s78.emit_stage(
        scripts,
        stage7_raw,
        current,
        out,
        "stage7_bounded_lookahead",
        args.segment_anchor_bases,
        timings,
    )

    stage8_paths, stage8_stats = s78.resolve_lookahead_seeded_paths(
        simplified,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        args.dominance,
        args.min_direct,
        args.min_length,
        args.lookahead_depth,
        args.lookahead_max_branch,
        args.lookahead_discount,
        args.lookahead_dominance,
        args.lookahead_margin,
    )
    stage8_raw = repeat_opt / "stage8_simplified_lookahead.fasta"
    stage8_write = gp.write_paths(
        stage8_paths,
        simplified,
        stage8_raw,
        repeat_opt / "stage8_simplified_lookahead.paths.tsv",
        args.min_length,
    )
    stage8 = s78.emit_stage(
        scripts,
        stage8_raw,
        current,
        out,
        "stage8_simplified_lookahead",
        args.segment_anchor_bases,
        timings,
    )

    backbone = out / "bridge_backbone.fasta"
    shutil.copy2(stage8, backbone)
    stats = {
        "pipeline": "bridge-backbone-v4-stage8",
        "promotion": {
            "production_backbone": "stage8_simplified_lookahead",
            "legacy_stage9_residual_unique": "disabled_not_promoted",
            "reason": "200k Zymo: negligible GF gain with misassemblies 11->30",
        },
        "lookahead": {
            "depth": args.lookahead_depth,
            "max_branch": args.lookahead_max_branch,
            "discount": args.lookahead_discount,
            "dominance": args.lookahead_dominance,
            "margin": args.lookahead_margin,
        },
        "full_read_threading": raw_stats,
        "projection_threading": projection_stats,
        "second_pass_threading": second_stats,
        "pair_repeat_threading": pair_stats,
        "graph_simplification": simplify_stats,
        "stage7_bounded_lookahead": {**stage7_stats, **stage7_write},
        "stage8_simplified_lookahead": {**stage8_stats, **stage8_write},
        "timings_seconds": timings,
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        / 1024.0,
        "outputs": {
            "stage7_bounded_lookahead": str(stage7),
            "stage8_simplified_lookahead": str(stage8),
            "production_backbone": str(backbone),
        },
    }
    (out / "stage8_backbone_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )

    manifest_path = out / "pipeline_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["pipeline"] = "bridge-backbone-v4-stage8"
    manifest["production_backbone"] = "stage8_simplified_lookahead"
    manifest["legacy_stage9"] = "disabled_not_promoted"
    manifest["stage8_backbone"] = stats
    outputs = manifest.setdefault("outputs", {})
    outputs.update(stats["outputs"])
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
