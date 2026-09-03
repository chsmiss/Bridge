#!/usr/bin/env python3
"""Stage 11 exploratory assembly: abundance-flow traversal + aggressive rare recall.

This stage is intentionally experimental and never replaces the promoted Stage8
backbone automatically. It combines ideas that have recently worked well in
metagenome assemblers:

* abundance-aware traversal: coverage flow may break otherwise unresolved graph
  ties, but only with independent physical/context evidence;
* deeper bounded lookahead over the existing graph (no invented graph edges);
* a wider residual-read multi-k ladder (k=15/17/21/25/31) with cross-k
  consensus, inspired by abundance-adaptive iterative assembly;
* pair-proposed positive gaps that are accepted only after sequence-resolved
  multi-k local assembly consensus. Unresolved N-gaps are split before output.

The script emits several Pareto candidates so the benchmark can separately
measure GF gain, continuity gain, and the cost in duplication/misassemblies.
"""
from __future__ import annotations

import argparse
import json
import math
import resource
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import graph_path_phaser as gp
import low_abundance_rescue as lr
import repeat_graph_optimizer as rg
import stage10_multik_rescue as s10
import stage789_optimizer as s78


@dataclass
class FlowRank:
    choice: gp.Choice
    score: float
    coverage_ratio: float
    sibling_share: float
    evidence_channels: int
    physical: int


def edge_evidence(graph: gp.Graph, history: list[str], candidate: str, forward: bool) -> gp.EdgeEvidence:
    edge = (history[-1], candidate) if forward else (candidate, history[0])
    return graph.edge.get(edge, gp.EdgeEvidence())


def evidence_channels(choice: gp.Choice, ev: gp.EdgeEvidence) -> int:
    channels = 0
    channels += int(ev.direct >= 2)
    channels += int(ev.gapped >= 2)
    channels += int(ev.pairs >= 2)
    channels += int(choice.raw_len >= 2 and choice.raw_support >= 2)
    channels += int(choice.second_len >= 3 and choice.second_support >= 1)
    channels += int(choice.proj_len >= 3 and choice.proj_support >= 1)
    channels += int(choice.high_len >= 3 and choice.high_support >= 1)
    return channels


def rank_flow_candidates(
    graph: gp.Graph,
    history: list[str],
    candidates: list[str],
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
) -> list[FlowRank]:
    available = [
        uid
        for uid in candidates
        if uid not in used and graph.rev.get(uid, uid) not in used
    ]
    if not available:
        return []
    current = history[-1] if forward else history[0]
    current_cov = max(0.001, graph.coverage.get(current, 0.0))
    total_cov = sum(max(0.001, graph.coverage.get(uid, 0.0)) for uid in available)
    ranked: list[FlowRank] = []
    for uid in available:
        choice = gp.candidate_choice(
            graph, history, uid, forward, raw_ctx, proj_ctx, high_ctx, repeat_ctx
        )
        ev = edge_evidence(graph, history, uid, forward)
        target_cov = max(0.001, graph.coverage.get(uid, 0.0))
        coverage_ratio = min(current_cov, target_cov) / max(current_cov, target_cov)
        sibling_share = target_cov / max(total_cov, 0.001)
        physical = max(ev.direct, ev.gapped, ev.pairs)
        channels = evidence_channels(choice, ev)
        strong = gp.strong_context(choice)

        # Flow is only a tie breaker. A coverage-compatible edge without read,
        # pair, projection, or high-k evidence is not allowed to create a path.
        if physical < 2 and not strong:
            continue
        if channels < 2:
            continue
        if coverage_ratio < 0.25:
            continue
        if sibling_share < 0.48 and (not strong or channels < 3):
            continue

        score = (
            max(0.0, choice.score)
            + 90.0 * coverage_ratio
            + 70.0 * min(1.0, sibling_share)
            + 18.0 * math.log1p(physical)
            + (35.0 if strong else 0.0)
            + 8.0 * min(channels, 5)
        )
        ranked.append(
            FlowRank(
                choice=choice,
                score=score,
                coverage_ratio=coverage_ratio,
                sibling_share=sibling_share,
                evidence_channels=channels,
                physical=physical,
            )
        )
    ranked.sort(
        key=lambda item: (
            -item.score,
            -item.evidence_channels,
            -item.coverage_ratio,
            -item.sibling_share,
            item.choice.uid,
        )
    )
    return ranked


def choose_extension_abundance_flow(
    graph: gp.Graph,
    history: list[str],
    candidates: list[str],
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
) -> tuple[gp.Choice | None, str]:
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
        0.66,
        3,
    )
    if direct is not None:
        return direct, "direct"

    lookahead, rescued = s78.choose_extension_lookahead(
        graph,
        history,
        candidates,
        used,
        forward,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        0.66,
        3,
        5,
        5,
        0.75,
        0.53,
        1.07,
    )
    if lookahead is not None:
        return lookahead, "lookahead" if rescued else "direct"

    ranked = rank_flow_candidates(
        graph,
        history,
        candidates,
        used,
        forward,
        raw_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
    )
    if not ranked:
        return None, "stop"
    best = ranked[0]
    second = ranked[1].score if len(ranked) > 1 else 0.0
    total = sum(max(0.0, item.score) for item in ranked)
    if total <= 0 or best.score / total < 0.52:
        return None, "stop"
    if second > 0 and best.score < second * 1.08:
        return None, "stop"
    return best.choice, "flow"


def resolve_abundance_flow_paths(
    graph: gp.Graph,
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
        if graph.ambiguous(seed):
            stats["ambiguous_paths_seeded_after_unique"] += 1
        path = [seed]
        local_seen = {seed, graph.rev.get(seed, seed)}
        while True:
            current = path[0]
            choice, mode = choose_extension_abundance_flow(
                graph,
                path,
                graph.inc.get(current, []),
                used | local_seen,
                False,
                raw_ctx,
                proj_ctx,
                high_ctx,
                repeat_ctx,
            )
            if choice is None:
                if len(graph.inc.get(current, [])) > 1:
                    stats["branch_stops"] += 1
                break
            path.insert(0, choice.uid)
            local_seen.update((choice.uid, graph.rev.get(choice.uid, choice.uid)))
            stats[f"{mode}_extensions"] += 1
        while True:
            current = path[-1]
            choice, mode = choose_extension_abundance_flow(
                graph,
                path,
                graph.out.get(current, []),
                used | local_seen,
                True,
                raw_ctx,
                proj_ctx,
                high_ctx,
                repeat_ctx,
            )
            if choice is None:
                if len(graph.out.get(current, [])) > 1:
                    stats["branch_stops"] += 1
                break
            path.append(choice.uid)
            local_seen.update((choice.uid, graph.rev.get(choice.uid, choice.uid)))
            stats[f"{mode}_extensions"] += 1
        gp.claim(path, graph, used)
        if len(gp.path_sequence(path, graph)) >= min_length:
            paths.append(path)
    stats["paths"] = len(paths)
    stats["claimed_orientations"] = len(used)
    return paths, dict(stats)


def run(cmd: list[object]) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True)
    return time.monotonic() - started


def positive_gap_localfill(
    scripts: Path,
    contigs: Path,
    read1: Path,
    read2: Path,
    outdir: Path,
    threads: int,
    timings: dict[str, float],
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    scaffold = outdir / "pair_scaffold.fasta"
    links = outdir / "pair_links.tsv"
    final = outdir / "primary_contigs.fasta"
    timings[f"pair_scaffold_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "pair_gap_refine.py",
            contigs,
            "-1",
            read1,
            "-2",
            read2,
            "-o",
            scaffold,
            "--links",
            links,
            "--threads",
            threads,
            "--min-mapq",
            30,
            "--min-support",
            4,
            "--dominance",
            0.92,
            "--end-window",
            650,
            "--min-overlap",
            31,
            "--max-gap",
            350,
        ]
    )
    timings[f"gap_consensus_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "fill_scaffold_gaps_multik.py",
            scaffold,
            "-1",
            read1,
            "-2",
            read2,
            "-o",
            final,
            "--report",
            outdir / "gap_consensus.tsv",
            "--anchor-k",
            31,
            "--local-ks",
            "17,21,25",
            "--min-consensus",
            2,
            "--dominance",
            0.68,
            "--flank",
            220,
            "--min-length",
            200,
        ]
    )
    return final


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--segment-anchor-bases", type=int, default=31)
    ap.add_argument("--ks", default="15,17,21,25,31")
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
    required = [
        current,
        backbone,
        target_gfa,
        projection_primary,
        highk_gfa,
        base_paths,
        args.read1,
        args.read2,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage11 inputs: " + ", ".join(missing))

    stage11 = out / "stage11_aggressive"
    stage11.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    # 1) Abundance-flow traversal over the simplified k31 graph.
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
    flow_paths, flow_stats = resolve_abundance_flow_paths(
        simplified, raw_ctx, proj_ctx, high_ctx, repeat_ctx, 200
    )
    flow_raw = stage11 / "flow_aggressive.raw.fasta"
    flow_write = gp.write_paths(
        flow_paths,
        simplified,
        flow_raw,
        stage11 / "flow_aggressive.paths.tsv",
        200,
    )
    flow_final = s78.emit_stage(
        scripts,
        flow_raw,
        current,
        out,
        "stage11_flow_aggressive",
        args.segment_anchor_bases,
        timings,
    )

    # 2) Wider low-abundance multi-k recall. k=15 is exploratory; candidate
    # promotion still requires agreement with other k values.
    ks = sorted({int(x) for x in args.ks.split(",") if x.strip()})
    if len(ks) < 3 or any(k < 15 or k > 63 for k in ks):
        raise SystemExit("--ks must contain >=3 k values in [15,63]")
    backbone31, backbone_bases = lr.backbone_kmers(backbone, 31)
    backbone21, _ = lr.backbone_kmers(backbone, 21)
    rare_r1 = stage11 / "rare_R1.fastq.gz"
    rare_r2 = stage11 / "rare_R2.fastq.gz"
    residual_pair_stats = lr.select_residual_pairs(
        args.read1,
        args.read2,
        rare_r1,
        rare_r2,
        backbone31,
        31,
        3,
        0.20,
        0.55,
    )
    inputs: dict[int, list[Path]] = {}
    if residual_pair_stats["pairs_kept"] >= 100:
        for k in ks:
            asm = stage11 / f"residual_k{k}"
            timings[f"rare_k{k}"] = s10.assemble_residual_k(
                args.bridgeasm, rare_r1, rare_r2, asm, k, args.threads
            )
            inputs[k] = [asm / "primary_contigs.fasta", asm / "haplotigs.fasta"]
    raw = s10.load_raw_candidates(inputs, 200)
    multik = s10.annotate_multik_candidates(raw, backbone31, backbone21)
    consensus = s10.select_multik_candidates(
        multik,
        backbone31,
        min_novel_kmers=48,
        min_novel_fraction=0.60,
        min_cross_sources=2,
        min_cross_fraction=0.30,
        max_total_bases=max(40_000, int(backbone_bases * 0.10)),
        max_fraction_per_k=0.55,
        allow_strong_single_k=False,
    )
    aggressive = s10.select_multik_candidates(
        multik,
        backbone31,
        min_novel_kmers=32,
        min_novel_fraction=0.45,
        min_cross_sources=1,
        min_cross_fraction=0.18,
        max_total_bases=max(60_000, int(backbone_bases * 0.15)),
        max_fraction_per_k=0.55,
        allow_strong_single_k=True,
    )
    consensus_add = stage11 / "rare_consensus_additions.fasta"
    aggressive_add = stage11 / "rare_aggressive_additions.fasta"
    lr.write_fasta(
        ((f"stage11_consensus_{i:06d}_k{c.k}", c.seq) for i, c in enumerate(consensus, 1)),
        consensus_add,
    )
    lr.write_fasta(
        ((f"stage11_aggressive_{i:06d}_k{c.k}", c.seq) for i, c in enumerate(aggressive, 1)),
        aggressive_add,
    )
    s10.write_metadata(consensus, stage11 / "rare_consensus.tsv")
    s10.write_metadata(aggressive, stage11 / "rare_aggressive.tsv")
    rare_consensus = lr.make_union_candidate(
        scripts,
        backbone,
        [consensus_add],
        stage11 / "candidate_rare_consensus",
        timings,
    )
    rare_aggressive = lr.make_union_candidate(
        scripts,
        backbone,
        [aggressive_add],
        stage11 / "candidate_rare_aggressive",
        timings,
    )
    flow_rare = lr.make_union_candidate(
        scripts,
        flow_final,
        [consensus_add],
        stage11 / "candidate_flow_rare_consensus",
        timings,
    )

    # 3) Allow positive pair gaps, but keep a join only if multi-k local
    # assembly resolves it. Any unresolved N is split back to contigs.
    stage8_local = positive_gap_localfill(
        scripts,
        backbone,
        args.read1,
        args.read2,
        stage11 / "candidate_stage8_localfill",
        args.threads,
        timings,
    )
    flow_local = positive_gap_localfill(
        scripts,
        flow_final,
        args.read1,
        args.read2,
        stage11 / "candidate_flow_localfill",
        args.threads,
        timings,
    )
    rare_local = positive_gap_localfill(
        scripts,
        rare_consensus,
        args.read1,
        args.read2,
        stage11 / "candidate_rare_consensus_localfill",
        args.threads,
        timings,
    )
    flow_rare_local = positive_gap_localfill(
        scripts,
        flow_rare,
        args.read1,
        args.read2,
        stage11 / "candidate_flow_rare_localfill",
        args.threads,
        timings,
    )

    outputs = {
        "stage8_backbone": str(backbone),
        "stage8_localfill": str(stage8_local),
        "flow_aggressive": str(flow_final),
        "flow_aggressive_localfill": str(flow_local),
        "rare_consensus": str(rare_consensus),
        "rare_consensus_localfill": str(rare_local),
        "rare_aggressive": str(rare_aggressive),
        "flow_rare_consensus": str(flow_rare),
        "flow_rare_consensus_localfill": str(flow_rare_local),
    }
    stats = {
        "pipeline": "bridge-stage11-aggressive-v1",
        "policy": {
            "production_backbone_replaced": False,
            "graph_edges_invented": False,
            "positive_pair_gap_requires_multik_sequence_consensus": True,
            "unresolved_n_gaps_split_before_output": True,
            "flow_override_min_evidence_channels": 2,
            "flow_lookahead_depth": 5,
            "residual_ks": ks,
        },
        "threading": {
            "raw": raw_stats,
            "projection": projection_stats,
            "second": second_stats,
            "pair": pair_stats,
        },
        "graph_simplification": simplify_stats,
        "flow": {**flow_stats, **flow_write},
        "residual_pairs": residual_pair_stats,
        "rare": {
            "raw_candidates": len(multik),
            "consensus_selected": len(consensus),
            "consensus_bases": sum(len(c.seq) for c in consensus),
            "aggressive_selected": len(aggressive),
            "aggressive_bases": sum(len(c.seq) for c in aggressive),
        },
        "timings_seconds": timings,
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0,
        "outputs": outputs,
    }
    (stage11 / "stage11_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    manifest_path = out / "pipeline_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["stage11_aggressive"] = stats
    manifest.setdefault("outputs", {}).update({f"stage11_{k}": v for k, v in outputs.items()})
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
