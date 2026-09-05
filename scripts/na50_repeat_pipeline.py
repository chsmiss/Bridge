#!/usr/bin/env python3
"""Cumulative NA50 pipeline with minimal-anchor segment emission and repeat optimization.

Stages 1-4 mirror graph_path_phaser.py. Stage 5 uses unique-flank-first repeat
traversal and stage 6 adds conservative graph simplification. Every stage passes
through the same boundary-aware segment replacement and short-graft postprocess,
so continuity gains are compared under a matched low-duplication policy.
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--segment-anchor-bases", type=int, default=31)
    args = ap.parse_args()
    if args.segment_anchor_bases < 31:
        raise SystemExit("segment-anchor-bases must be >=31")

    scripts = Path(__file__).resolve().parent
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}
    started = time.monotonic()

    base = out / "current_pipeline"
    timings["current_pipeline"] = run([
        sys.executable,
        scripts / "bridge_evidence_pipeline.py",
        "--bridgeasm", args.bridgeasm,
        "--read1", args.read1,
        "--read2", args.read2,
        "--output", base,
        "--threads", args.threads,
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
        "--gfa", target_gfa,
        "-1", args.read1,
        "-2", args.read2,
        "--projection", projection_primary,
        "--projection", projection_haplotigs,
        "--highk-gfa", highk_gfa,
        "-o", graph_opt,
        "--anchor-k", 31,
        "--max-context", 6,
        "--dominance", 0.72,
        "--min-direct", 4,
        "--min-length", 200,
    ])

    repeat_opt = out / "repeat_optimizer"
    timings["repeat_optimizer"] = run([
        sys.executable,
        scripts / "repeat_graph_optimizer_v2.py",
        "--gfa", target_gfa,
        "-1", args.read1,
        "-2", args.read2,
        "--base-paths", graph_opt / "stage4_second_pass.paths.tsv",
        "--projection", projection_primary,
        "--projection", projection_haplotigs,
        "--highk-gfa", highk_gfa,
        "-o", repeat_opt,
        "--anchor-k", 31,
        "--max-context", 8,
        "--max-pair-bridge-edges", 6,
        "--max-pair-span", 320,
        "--dominance", 0.70,
        "--min-direct", 4,
        "--min-length", 200,
    ])

    candidates = {
        "stage1_full_read": graph_opt / "stage1_full_read.fasta",
        "stage2_iterative_projection": graph_opt / "stage2_iterative_projection.fasta",
        "stage3_local_highk": graph_opt / "stage3_local_highk.fasta",
        "stage4_second_pass": graph_opt / "stage4_second_pass.fasta",
        "stage5_repeat_traversal": repeat_opt / "stage5_repeat_seeded.fasta",
        "stage6_graph_simplified": repeat_opt / "stage6_graph_simplified_seeded.fasta",
    }
    outputs: dict[str, str] = {"current": str(current)}
    for name, candidate in candidates.items():
        outputs[name] = str(
            emit_stage(
                scripts,
                candidate,
                current,
                out,
                name,
                args.segment_anchor_bases,
                timings,
            )
        )

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest = {
        "pipeline": "bridge-na50-repeat-v2-minimal-anchor-unique-flank",
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": usage.ru_maxrss / 1024.0,
        "segment_anchor_bases": args.segment_anchor_bases,
        "short_min_length": 31,
        "repeat_optimizer": "v2_unique_flank_seeded",
        "timings_seconds": timings,
        "outputs": outputs,
    }
    (out / "pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
