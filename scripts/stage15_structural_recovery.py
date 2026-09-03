#!/usr/bin/env python3
"""Stage 15: structural recovery for Bridge's measured NA50/GF bottlenecks.

1. Soft graph threading: re-thread the full library against the same k31
   unitig graph with unique k17/k21 anchors.  Only path-context support absent
   from exact k31 threading is added to the conservative Stage8 resolver.
2. Rare-sequence reintegration: Stage10 strict additions are cross-k validated
   but mostly below the old 500 bp virtualization cutoff.  Project those short
   trusted sequences back before k31 graph construction, then promote rebuilt
   primary paths through k41 and k55.
3. Re-run long-component abundance flow on the reintegrated k31 graph.  Final
   sequence joining remains exact (>=81 bp) and reference-free.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import graph_path_phaser as gp
import low_abundance_rescue as lr
import repeat_graph_optimizer as rg
import stage13_three_methods as s13
import stage14_amplified_methods as s14
import stage789_optimizer as s78


def run(cmd: list[object], *, env: dict[str, str] | None = None) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True, env=env)
    return time.monotonic() - started


def concat_gzip(inputs: list[Path], output: Path) -> None:
    """Concatenate gzip members without recompressing."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as out:
        for path in inputs:
            if not path.exists():
                continue
            with path.open("rb") as handle:
                shutil.copyfileobj(handle, out, 1 << 20)


def soft_context_increment(
    baseline: Counter[tuple[str, ...]],
    extras: list[Counter[tuple[str, ...]]],
    *,
    max_weight: int = 4,
) -> Counter[tuple[str, ...]]:
    """Return conservative lower-k context support absent from k31 threading.

    The same fragment can be rediscovered at k17 and k21, so take the maximum
    excess support across k values instead of summing it.  Two-node contexts
    need two extra observations; paths spanning >=2 edges need one because the
    linkage itself is already more specific.
    """
    result: Counter[tuple[str, ...]] = Counter()
    for counter in extras:
        for key, support in counter.items():
            if len(key) < 2:
                continue
            delta = max(0, support - baseline.get(key, 0))
            minimum = 1 if len(key) >= 3 else 2
            if delta < minimum:
                continue
            result[key] = max(result.get(key, 0), min(max_weight, delta))
    return result


def build_soft_thread_candidate(
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
    projection_primary = base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = base / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"
    base_paths = graph_opt / "stage4_second_pass.paths.tsv"
    required = [target_gfa, projection_primary, highk_gfa, base_paths, strict_baseline]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing soft-thread inputs: " + ", ".join(missing))

    graph = gp.Graph.from_gfa(target_gfa)
    membership = gp.preliminary_membership(rg.load_paths(base_paths))
    indexes = {k: gp.KmerIndex(graph, k) for k in (17, 21, 31)}

    ctx_by_k: dict[int, Counter[tuple[str, ...]]] = {}
    stats_by_k: dict[int, dict[str, int]] = {}
    started = time.monotonic()
    for k in (31, 21, 17):
        ctx, stats = gp.collect_read_contexts(graph, indexes[k], read1, read2, None, 10)
        ctx_by_k[k] = ctx
        stats_by_k[k] = stats
    increment = soft_context_increment(ctx_by_k[31], [ctx_by_k[21], ctx_by_k[17]])
    enhanced_raw = Counter(ctx_by_k[31])
    enhanced_raw.update(increment)
    timings["soft_thread_full_library"] = time.monotonic() - started

    proj_ctx, high_ctx, projection_stats = rg.collect_projection_contexts(
        graph,
        indexes[31],
        [projection_primary, projection_haplotigs],
        [highk_gfa],
        repeat_opt,
        8,
    )
    second_ctx, second_stats = gp.collect_read_contexts(
        graph, indexes[31], read1, read2, membership, 8
    )
    for key in list(second_ctx):
        baseline = ctx_by_k[31].get(key, 0)
        if second_ctx[key] <= baseline:
            del second_ctx[key]
        else:
            second_ctx[key] -= baseline
    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph, indexes[31], read1, read2, membership, 8, 8, 420
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
    stage15_dir = pipeline_dir / "stage15_structural_recovery" / "soft_thread"
    stage15_dir.mkdir(parents=True, exist_ok=True)
    raw_fasta = stage15_dir / "soft_thread_raw.fasta"
    write_stats = gp.write_paths(
        paths,
        simplified,
        raw_fasta,
        stage15_dir / "soft_thread.paths.tsv",
        200,
    )
    final = s78.emit_stage(
        scripts,
        raw_fasta,
        strict_baseline,
        pipeline_dir,
        "stage15_soft_thread",
        31,
        timings,
    )
    copied = stage15_dir / "candidate_soft_thread" / "primary_contigs.fasta"
    copied.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(final, copied)

    base_threaded = stats_by_k[31]["threaded_reads"]
    soft_threaded = max(stats_by_k[17]["threaded_reads"], stats_by_k[21]["threaded_reads"])
    return copied, {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "k31_threading": stats_by_k[31],
        "k21_threading": stats_by_k[21],
        "k17_threading": stats_by_k[17],
        "threaded_read_gain_vs_k31": soft_threaded - base_threaded,
        "incremental_contexts": len(increment),
        "incremental_context_weight": sum(increment.values()),
        "projection": projection_stats,
        "second_pass": second_stats,
        "pair_threading": pair_stats,
        "simplification": simplify_stats,
        "path_resolution": resolve_stats,
        "write": write_stats,
    }


def virtualize(
    scripts: Path,
    fasta_inputs: list[Path],
    out1: Path,
    out2: Path,
    *,
    read_length: int,
    insert_size: int,
    stride: int,
    min_length: int,
) -> float:
    return run(
        [
            sys.executable,
            scripts / "make_virtual_pairs.py",
            *fasta_inputs,
            "--read1",
            out1,
            "--read2",
            out2,
            "--read-length",
            read_length,
            "--insert-size",
            insert_size,
            "--stride",
            stride,
            "--min-length",
            min_length,
            "--max-pairs-per-record",
            32,
        ]
    )


def run_reintegration_assembly(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    outdir: Path,
    k: int,
    threads: int,
    *,
    protect_sparse: bool,
) -> float:
    mercy = {31: 24, 41: 16, 55: 12}.get(k, 12)
    env = os.environ.copy()
    if protect_sparse:
        env["BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION"] = "0.85"
        env["BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY"] = "37"
        env["BRIDGEASM_MATE_TERMINAL_MERCY_KMERS"] = "64"
    else:
        env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION", None)
        env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY", None)
        env.pop("BRIDGEASM_MATE_TERMINAL_MERCY_KMERS", None)
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
            4,
            "--primary-dominance",
            0.78,
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


def profile(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.exists() else {}


def build_rare_reintegration_candidate(
    scripts: Path,
    bridgeasm: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    read1: Path,
    read2: Path,
    threads: int,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    stage10 = pipeline_dir / "stage10_multik_rescue"
    strict_add = stage10 / "multik_strict_additions.fasta"
    outdir = pipeline_dir / "stage15_structural_recovery" / "rare_reintegration"
    outdir.mkdir(parents=True, exist_ok=True)

    rare_v1 = outdir / "rare.virtual_R1.fastq.gz"
    rare_v2 = outdir / "rare.virtual_R2.fastq.gz"
    timings["reintegrate_virtualize_rare"] = virtualize(
        scripts,
        [strict_add],
        rare_v1,
        rare_v2,
        read_length=91,
        insert_size=190,
        stride=60,
        min_length=190,
    )
    aug31_1, aug31_2 = outdir / "k31.aug_R1.fastq.gz", outdir / "k31.aug_R2.fastq.gz"
    concat_gzip([read1, rare_v1], aug31_1)
    concat_gzip([read2, rare_v2], aug31_2)
    k31 = outdir / "k31_reintegrated"
    timings["reintegrate_k31"] = run_reintegration_assembly(
        bridgeasm, aug31_1, aug31_2, k31, 31, threads, protect_sparse=True
    )

    v31_1, v31_2 = outdir / "k31.virtual_R1.fastq.gz", outdir / "k31.virtual_R2.fastq.gz"
    timings["reintegrate_virtualize_k31"] = virtualize(
        scripts,
        [k31 / "primary_contigs.fasta", k31 / "haplotigs.fasta"],
        v31_1,
        v31_2,
        read_length=101,
        insert_size=250,
        stride=120,
        min_length=250,
    )
    aug41_1, aug41_2 = outdir / "k41.aug_R1.fastq.gz", outdir / "k41.aug_R2.fastq.gz"
    concat_gzip([read1, rare_v1, v31_1], aug41_1)
    concat_gzip([read2, rare_v2, v31_2], aug41_2)
    k41 = outdir / "k41_promoted"
    timings["reintegrate_k41"] = run_reintegration_assembly(
        bridgeasm, aug41_1, aug41_2, k41, 41, threads, protect_sparse=False
    )

    v41_1, v41_2 = outdir / "k41.virtual_R1.fastq.gz", outdir / "k41.virtual_R2.fastq.gz"
    timings["reintegrate_virtualize_k41"] = virtualize(
        scripts,
        [k41 / "primary_contigs.fasta", k41 / "haplotigs.fasta"],
        v41_1,
        v41_2,
        read_length=101,
        insert_size=250,
        stride=120,
        min_length=250,
    )
    aug55_1, aug55_2 = outdir / "k55.aug_R1.fastq.gz", outdir / "k55.aug_R2.fastq.gz"
    concat_gzip([read1, rare_v1, v41_1], aug55_1)
    concat_gzip([read2, rare_v2, v41_2], aug55_2)
    k55 = outdir / "k55_promoted"
    timings["reintegrate_k55"] = run_reintegration_assembly(
        bridgeasm, aug55_1, aug55_2, k55, 55, threads, protect_sparse=False
    )

    additions = [
        k31 / "primary_contigs.fasta",
        k41 / "primary_contigs.fasta",
        k55 / "primary_contigs.fasta",
    ]
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        additions,
        outdir / "candidate_reintegrated",
        timings,
        min_overlap=81,
    )
    return final, k31 / "assembly.gfa", {
        "rare_source_records": sum(1 for _ in lr.fasta_records(strict_add)),
        "k31": profile(k31 / "run_profile.json"),
        "k41": profile(k41 / "run_profile.json"),
        "k55": profile(k55 / "run_profile.json"),
        "policy": {
            "rare_virtual_min_length": 190,
            "k31_sparse_protection": "singleton>=0.85 Q>=37 plus mate-terminal mercy 64",
            "higher_k": "raw + trusted rare virtual + previous-stage virtual paths",
            "final_join": "exact overlap >=81 bp",
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--ks", default="17,21,25,31")
    args = ap.parse_args()

    started = time.monotonic()
    scripts = Path(__file__).resolve().parent
    out = args.pipeline_dir
    stage10 = out / "stage10_multik_rescue"
    strict_baseline = stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    backbone = out / "bridge_backbone.fasta"
    required = [
        args.bridgeasm,
        args.read1,
        args.read2,
        strict_baseline,
        backbone,
        stage10 / "multik_strict_additions.fasta",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage15 inputs: " + ", ".join(missing))

    stage15 = out / "stage15_structural_recovery"
    stage15.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    soft_final, soft_stats = build_soft_thread_candidate(
        scripts, out, strict_baseline, args.read1, args.read2, timings
    )
    reintegrated_final, reintegrated_gfa, reintegrated_stats = build_rare_reintegration_candidate(
        scripts,
        args.bridgeasm,
        out,
        strict_baseline,
        args.read1,
        args.read2,
        args.threads,
        timings,
    )

    ks = sorted({int(value) for value in args.ks.split(",") if value.strip()})
    inputs = s13.residual_inputs(stage10, ks)
    pools = s13.cross_k_pools(inputs)
    flow_final, flow_add, flow_stats = s14.build_long_flow_candidate(
        scripts,
        strict_baseline,
        reintegrated_gfa,
        args.read1,
        args.read2,
        pools,
        stage15 / "reintegrated_flow",
        timings,
    )
    reintegrated_flow = s14.make_bridge_candidate(
        scripts,
        reintegrated_final,
        [flow_add],
        stage15 / "reintegrated_plus_flow" / "candidate_reintegrated_plus_flow",
        timings,
        min_overlap=81,
    )
    combined = s14.make_bridge_candidate(
        scripts,
        soft_final,
        [reintegrated_final, flow_add],
        stage15 / "combined" / "candidate_combined",
        timings,
        min_overlap=81,
    )

    stats = {
        "pipeline": "bridge-stage15-structural-recovery-v1",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "soft_threading": "unique k17/k21 anchors on k31 unitig graph; max excess context only",
            "rare_reintegration": "Stage10 strict additions virtualized before k31 graph build",
            "persistent_multik": "k31 -> virtual paths -> k41 -> virtual paths -> k55",
            "long_flow": "long-component NNLS rerun on reintegrated k31 graph",
            "sequence_join": "exact overlap >=81 bp",
        },
        "methods": {
            "soft_thread": soft_stats,
            "rare_reintegration": reintegrated_stats,
            "reintegrated_flow": flow_stats,
        },
        "outputs": {
            "soft_thread": str(soft_final),
            "rare_reintegration": str(reintegrated_final),
            "flow_on_reintegrated_graph": str(flow_final),
            "reintegrated_plus_flow": str(reintegrated_flow),
            "combined": str(combined),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (stage15 / "stage15_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
