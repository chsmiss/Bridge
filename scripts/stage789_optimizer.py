#!/usr/bin/env python3
"""Append conservative stages 7-9 to an existing NA50 repeat pipeline run.

Stage 7 adds bounded lookahead only at ambiguous junctions rejected by the
existing one-step chooser. Stage 8 repeats that traversal after conservative
graph simplification. Stage 9 adds only strongly supported residual sequence
patches carrying novel 31-mers. No stage invents graph edges; stage 7/8 paths
claim each oriented node at most once, and all three candidates pass through
the same low-duplication backbone replacement and graft postprocess used by the
stage1-6 pipeline.
"""
from __future__ import annotations

import argparse
import json
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import graph_path_phaser as gp
import repeat_graph_optimizer as rg


def run(cmd) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True)
    return time.monotonic() - started


def emit_stage(
    scripts: Path,
    candidate: Path,
    recovery: Path,
    out: Path,
    name: str,
    anchor_bases: int,
    timings: dict[str, float],
) -> Path:
    stage_dir = out / name
    stage_dir.mkdir(parents=True, exist_ok=True)
    replacement = stage_dir / "backbone_replacement_union.fasta"
    timings[f"replace_{name}"] = run([
        sys.executable,
        scripts / "merge_backbone_replacement.py",
        "--backbone", candidate,
        "--recovery", recovery,
        "-o", replacement,
        "--report", stage_dir / "backbone_replacement.tsv",
        "--stats-json", stage_dir / "backbone_replacement.json",
        "-k", 31,
        "--replace-fraction", 0.85,
        "--min-informative-kmers", 20,
        "--segment-anchor-bases", anchor_bases,
        "--min-novel-kmers", 4,
        "--merge-represented-gap-kmers", 8,
        "--max-novel-hole-kmers", 2,
        "--min-segment-length", 31,
    ])
    graft = stage_dir / "graft"
    timings[f"graft_{name}"] = run([
        sys.executable,
        scripts / "postprocess_segment_grafts.py",
        replacement,
        "--output-dir", graft,
        "--short-min-length", 31,
        "--final-min-length", 200,
        "--min-overlap", 31,
        "--overlap-margin", 10,
    ])
    final = graft / "primary_contigs.fasta"
    copied = out / f"{name}.fasta"
    shutil.copy2(final, copied)
    return copied


def concatenate_fastas(inputs: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        for path in inputs:
            if not path.exists() or path.stat().st_size == 0:
                continue
            data = path.read_bytes()
            handle.write(data)
            if data and not data.endswith(b"\n"):
                handle.write(b"\n")


@dataclass
class LookaheadBranch:
    choice: gp.Choice
    score: float
    context_hits: int
    steps: int


def advance_history(history: list[str], candidate: str, forward: bool) -> list[str]:
    return history + [candidate] if forward else [candidate] + history


def next_candidates(graph: gp.Graph, uid: str, forward: bool) -> list[str]:
    return graph.out.get(uid, []) if forward else graph.inc.get(uid, [])


def physical_support(
    graph: gp.Graph, history: list[str], candidate: str, forward: bool
) -> int:
    edge = (history[-1], candidate) if forward else (candidate, history[0])
    ev = graph.edge.get(edge, gp.EdgeEvidence())
    return max(ev.direct, ev.gapped, ev.pairs)


def lookahead_tail_score(
    graph: gp.Graph,
    history: list[str],
    current: str,
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
    depth: int,
    max_branch: int,
    discount: float,
) -> tuple[float, int, int]:
    if depth <= 0:
        return 0.0, 0, 0
    available = [
        uid
        for uid in next_candidates(graph, current, forward)
        if uid not in used and graph.rev.get(uid, uid) not in used
    ]
    if not available:
        return 0.0, 0, 0
    ranked = []
    for uid in available:
        choice = gp.candidate_choice(
            graph, history, uid, forward, raw_ctx, proj_ctx, high_ctx, repeat_ctx
        )
        ranked.append((max(0.0, choice.score), uid, choice))
    ranked.sort(key=lambda item: (-item[0], item[1]))

    best = (0.0, 0, 0)
    for immediate, uid, choice in ranked[:max_branch]:
        history2 = advance_history(history, uid, forward)
        used2 = used | {uid, graph.rev.get(uid, uid)}
        tail_score, tail_hits, tail_steps = lookahead_tail_score(
            graph,
            history2,
            uid,
            used2,
            forward,
            raw_ctx,
            proj_ctx,
            high_ctx,
            repeat_ctx,
            depth - 1,
            max_branch,
            discount,
        )
        candidate = (
            immediate + discount * tail_score,
            int(gp.strong_context(choice)) + tail_hits,
            1 + tail_steps,
        )
        if candidate > best:
            best = candidate
    return best


def score_lookahead_branch(
    graph: gp.Graph,
    history: list[str],
    candidate: str,
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
    depth: int,
    max_branch: int,
    discount: float,
) -> LookaheadBranch:
    choice = gp.candidate_choice(
        graph, history, candidate, forward, raw_ctx, proj_ctx, high_ctx, repeat_ctx
    )
    history2 = advance_history(history, candidate, forward)
    used2 = used | {candidate, graph.rev.get(candidate, candidate)}
    tail_score, tail_hits, tail_steps = lookahead_tail_score(
        graph,
        history2,
        candidate,
        used2,
        forward,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        max(0, depth - 1),
        max_branch,
        discount,
    )
    return LookaheadBranch(
        choice=choice,
        score=max(0.0, choice.score) + discount * tail_score,
        context_hits=int(gp.strong_context(choice)) + tail_hits,
        steps=1 + tail_steps,
    )


def choose_extension_lookahead(
    graph: gp.Graph,
    history: list[str],
    candidates: list[str],
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
    dominance: float,
    min_direct: int,
    lookahead_depth: int,
    lookahead_max_branch: int,
    lookahead_discount: float,
    lookahead_dominance: float,
    lookahead_margin: float,
) -> tuple[gp.Choice | None, bool]:
    direct = gp.choose_extension(
        graph,
        history,
        candidates,
        used,
        forward,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        dominance,
        min_direct,
    )
    if direct is not None:
        return direct, False

    available = [
        uid
        for uid in candidates
        if uid not in used and graph.rev.get(uid, uid) not in used
    ]
    if len(available) <= 1 or lookahead_depth <= 1:
        return None, False

    branches = [
        score_lookahead_branch(
            graph,
            history,
            uid,
            used,
            forward,
            raw_ctx,
            proj_ctx,
            high_ctx,
            repeat_ctx,
            lookahead_depth,
            lookahead_max_branch,
            lookahead_discount,
        )
        for uid in available
    ]
    branches.sort(
        key=lambda item: (-item.score, -item.context_hits, -item.steps, item.choice.uid)
    )
    if not branches or branches[0].score <= 0:
        return None, False

    best = branches[0]
    second_score = branches[1].score if len(branches) > 1 else 0.0
    total = sum(max(0.0, item.score) for item in branches)
    fraction = best.score / total if total > 0 else 0.0
    if fraction < lookahead_dominance:
        return None, False
    if second_score > 0 and best.score < second_score * lookahead_margin:
        return None, False

    root_support = physical_support(graph, history, best.choice.uid, forward)
    if root_support < max(1, min_direct // 2) and not gp.strong_context(best.choice):
        return None, False
    return best.choice, True


def resolve_lookahead_seeded_paths(
    graph: gp.Graph,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
    dominance: float,
    min_direct: int,
    min_length: int,
    lookahead_depth: int,
    lookahead_max_branch: int,
    lookahead_discount: float,
    lookahead_dominance: float,
    lookahead_margin: float,
) -> tuple[list[list[str]], dict[str, int]]:
    used: set[str] = set()
    paths: list[list[str]] = []
    phased_extensions = 0
    lookahead_extensions = 0
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
    ambiguous_seeds = 0
    for seed in seeds:
        if seed in used or graph.rev.get(seed, seed) in used:
            continue
        if graph.ambiguous(seed):
            ambiguous_seeds += 1
        path = [seed]
        local_seen = {seed, graph.rev.get(seed, seed)}
        while True:
            current = path[0]
            choice, rescued = choose_extension_lookahead(
                graph,
                path,
                graph.inc.get(current, []),
                used | local_seen,
                False,
                raw_ctx,
                proj_ctx,
                high_ctx,
                repeat_ctx,
                dominance,
                min_direct,
                lookahead_depth,
                lookahead_max_branch,
                lookahead_discount,
                lookahead_dominance,
                lookahead_margin,
            )
            if choice is None:
                if len(graph.inc.get(current, [])) > 1:
                    branch_stops += 1
                break
            path.insert(0, choice.uid)
            local_seen.update((choice.uid, graph.rev.get(choice.uid, choice.uid)))
            lookahead_extensions += int(rescued)
            phased_extensions += int(gp.strong_context(choice))
        while True:
            current = path[-1]
            choice, rescued = choose_extension_lookahead(
                graph,
                path,
                graph.out.get(current, []),
                used | local_seen,
                True,
                raw_ctx,
                proj_ctx,
                high_ctx,
                repeat_ctx,
                dominance,
                min_direct,
                lookahead_depth,
                lookahead_max_branch,
                lookahead_discount,
                lookahead_dominance,
                lookahead_margin,
            )
            if choice is None:
                if len(graph.out.get(current, [])) > 1:
                    branch_stops += 1
                break
            path.append(choice.uid)
            local_seen.update((choice.uid, graph.rev.get(choice.uid, choice.uid)))
            lookahead_extensions += int(rescued)
            phased_extensions += int(gp.strong_context(choice))
        gp.claim(path, graph, used)
        if len(gp.path_sequence(path, graph)) >= min_length:
            paths.append(path)
    return paths, {
        "paths": len(paths),
        "phased_extensions": phased_extensions,
        "lookahead_extensions": lookahead_extensions,
        "branch_stops": branch_stops,
        "claimed_orientations": len(used),
        "ambiguous_paths_seeded_after_unique": ambiguous_seeds,
    }


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
    all_ctx = rg.combined_contexts(raw_ctx, proj_ctx, high_ctx, repeat_ctx)
    simplified, simplify_stats = rg.simplify_graph(graph, all_ctx)

    stage7_paths, stage7_stats = resolve_lookahead_seeded_paths(
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
    stage7 = emit_stage(
        scripts, stage7_raw, current, out, "stage7_bounded_lookahead",
        args.segment_anchor_bases, timings
    )

    stage8_paths, stage8_stats = resolve_lookahead_seeded_paths(
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
    stage8 = emit_stage(
        scripts, stage8_raw, current, out, "stage8_simplified_lookahead",
        args.segment_anchor_bases, timings
    )

    residual_dir = out / "residual_rescue"
    residual_dir.mkdir(parents=True, exist_ok=True)
    residual = residual_dir / "stage9_residual_patches.fasta"
    residual_meta = residual_dir / "stage9_residual_patches.tsv"
    timings["stage9_residual_extract"] = run([
        sys.executable,
        scripts / "residual_path_cover.py",
        target_gfa,
        highk_gfa,
        "--backbone", stage8,
        "-o", residual,
        "--metadata", residual_meta,
        "--secondary-dominance", 0.25,
        "--extension-dominance", 0.80,
        "--min-support", 4,
        "--max-copy", 2,
        "--novel-k", 31,
        "--flank", 80,
        "--max-novel-gap", 64,
        "--min-novel-kmers", 8,
        "--min-novel-fraction", 0.10,
        "--min-length", 200,
        "--max-patch-length", 1000,
        "--max-patches", 500,
        "--max-total-fraction", 0.05,
    ])
    stage9_candidate = residual_dir / "stage9_residual_candidate.fasta"
    concatenate_fastas([stage8, residual], stage9_candidate)
    stage9 = emit_stage(
        scripts, stage9_candidate, current, out, "stage9_residual_unique",
        args.segment_anchor_bases, timings
    )

    stats = {
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
        "stage9_residual_unique": {
            "max_total_fraction": 0.05,
            "min_novel_fraction": 0.10,
            "min_novel_kmers": 8,
            "min_support": 4,
            "max_copy": 2,
        },
        "timings_seconds": timings,
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0,
        "outputs": {
            "stage7_bounded_lookahead": str(stage7),
            "stage8_simplified_lookahead": str(stage8),
            "stage9_residual_unique": str(stage9),
        },
    }
    (out / "stage789_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )

    manifest_path = out / "pipeline_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        manifest["pipeline"] = "bridge-na50-repeat-v3-stage789"
        manifest["repeat_optimizer"] = "v3_bounded_lookahead_residual"
        manifest["stage789"] = stats
        outputs = manifest.setdefault("outputs", {})
        outputs.update(stats["outputs"])
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
