#!/usr/bin/env python3
"""Cumulative NA50 pipeline with segment emission and repeat optimization.

Stages 1-4 mirror graph_path_phaser.py. Stage 5 adds mate-pair repeat traversal
and stage 6 adds conservative graph simplification. Every stage is emitted via
the same segment-level backbone replacement so continuity gains are compared at
matched duplication policy.
"""
from __future__ import annotations

import argparse
import json
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path


def run(cmd) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True)
    return time.monotonic() - started


def postprocess(scripts: Path, input_fasta: Path, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    merged = outdir / "union.fasta"
    filtered = outdir / "noncontained.fasta"
    final = outdir / "primary_contigs.fasta"
    run([
        sys.executable,
        scripts / "merge_fasta_unique.py",
        merged,
        input_fasta,
        "--min-length",
        200,
    ])
    run([
        sys.executable,
        scripts / "filter_contained_fasta.py",
        merged,
        filtered,
        "--min-length",
        200,
        "--seed-k",
        21,
        "--window",
        12,
        "--candidate-minimizers",
        16,
        "--removed-tsv",
        outdir / "contained_removed.tsv",
        "--stats-json",
        outdir / "containment_stats.json",
    ])
    run([
        sys.executable,
        scripts / "stitch_exact_overlaps.py",
        final,
        filtered,
        "--min-overlap",
        31,
        "--overlap-margin",
        10,
        "--seed-length",
        31,
        "--max-seed-occurrences",
        64,
        "--min-length",
        200,
    ])
    return final


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
    t0 = time.monotonic()
    timings[f"replace_{name}"] = run([
        sys.executable,
        scripts / "merge_backbone_replacement.py",
        "--backbone",
        candidate,
        "--recovery",
        recovery,
        "-o",
        replacement,
        "--report",
        stage_dir / "backbone_replacement.tsv",
        "--stats-json",
        stage_dir / "backbone_replacement.json",
        "-k",
        31,
        "--replace-fraction",
        0.85,
        "--min-informative-kmers",
        20,
        "--segment-anchor-bases",
        anchor_bases,
        "--min-novel-kmers",
        4,
        "--merge-represented-gap-kmers",
        8,
        "--max-novel-hole-kmers",
        2,
        "--min-segment-length",
        200,
    ])
    final = postprocess(scripts, replacement, stage_dir)
    timings[f"postprocess_{name}"] = time.monotonic() - t0
    copied = out / f"{name}.fasta"
    shutil.copy2(final, copied)
    return copied


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--segment-anchor-bases", type=int, default=64)
    args = ap.parse_args()

    scripts = Path(__file__).resolve().parent
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}
    started = time.monotonic()

    base = out / "current_pipeline"
    timings["current_pipeline"] = run([
        sys.executable,
        scripts / "bridge_evidence_pipeline.py",
        "--bridgeasm",
        args.bridgeasm,
        "--read1",
        args.read1,
        "--read2",
        args.read2,
        "--output",
        base,
        "--threads",
        args.threads,
    ])
    current = base / "step6_strain_projection.fasta"
    target_gfa = base / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = base / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"

    graph_opt = out / "graph_optimizer"
    timings["graph_optimizer"] = run([
        sys.executable,
        scripts / "graph_path_phaser.py",
        "--gfa",
        target_gfa,
        "-1",
        args.read1,
        "-2",
        args.read2,
        "--projection",
        projection_primary,
        "--projection",
        projection_haplotigs,
        "--highk-gfa",
        highk_gfa,
        "-o",
        graph_opt,
        "--anchor-k",
        31,
        "--max-context",
        6,
        "--dominance",
        0.72,
        "--min-direct",
        4,
        "--min-length",
        200,
    ])

    repeat_opt = out / "repeat_optimizer"
    timings["repeat_optimizer"] = run([
        sys.executable,
        scripts / "repeat_graph_optimizer.py",
        "--gfa",
        target_gfa,
        "-1",
        args.read1,
        "-2",
        args.read2,
        "--base-paths",
        graph_opt / "stage4_second_pass.paths.tsv",
        "--projection",
        projection_primary,
        "--projection",
        projection_haplotigs,
        "--highk-gfa",
        highk_gfa,
        "-o",
        repeat_opt,
        "--anchor-k",
        31,
        "--max-context",
        8,
        "--max-pair-bridge-edges",
        6,
        "--max-pair-span",
        320,
        "--dominance",
        0.70,
        "--min-direct",
        4,
        "--min-length",
        200,
    ])

    candidates = {
        "stage1_full_read": graph_opt / "stage1_full_read.fasta",
        "stage2_iterative_projection": graph_opt / "stage2_iterative_projection.fasta",
        "stage3_local_highk": graph_opt / "stage3_local_highk.fasta",
        "stage4_second_pass": graph_opt / "stage4_second_pass.fasta",
        "stage5_repeat_traversal": repeat_opt / "stage5_repeat_traversal.fasta",
        "stage6_graph_simplified": repeat_opt / "stage6_graph_simplified.fasta",
    }
    outputs: dict[str, str] = {"current": str(current)}
    for name, candidate in candidates.items():
        final = emit_stage(
            scripts,
            candidate,
            current,
            out,
            name,
            args.segment_anchor_bases,
            timings,
        )
        outputs[name] = str(final)

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest = {
        "pipeline": "bridge-na50-repeat-v1",
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": usage.ru_maxrss / 1024.0,
        "segment_anchor_bases": args.segment_anchor_bases,
        "timings_seconds": timings,
        "outputs": outputs,
    }
    (out / "pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
