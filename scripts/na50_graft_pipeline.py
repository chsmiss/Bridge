#!/usr/bin/env python3
"""NA50 graph pipeline with short-segment graft-before-filter emission.

This keeps the proven graph path stages unchanged. Recovery segments may be as
short as 64 bp and are allowed to exact-stitch into the graph backbone before
the final 200 bp contig cutoff, avoiding represented-sequence padding solely to
survive early filtering.
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


def emit(
    scripts: Path,
    candidate: Path,
    recovery: Path,
    out: Path,
    name: str,
    timings: dict[str, float],
) -> Path:
    stage = out / name
    stage.mkdir(parents=True, exist_ok=True)
    replacement = stage / "backbone_replacement_union.fasta"
    timings[f"replace_{name}"] = run([
        sys.executable,
        scripts / "merge_backbone_replacement.py",
        "--backbone", candidate,
        "--recovery", recovery,
        "-o", replacement,
        "--report", stage / "backbone_replacement.tsv",
        "--stats-json", stage / "backbone_replacement.json",
        "-k", 31,
        "--replace-fraction", 0.85,
        "--min-informative-kmers", 20,
        "--segment-anchor-bases", 64,
        "--min-novel-kmers", 4,
        "--merge-represented-gap-kmers", 8,
        "--max-novel-hole-kmers", 2,
        "--min-segment-length", 64,
    ])
    graft = stage / "graft"
    timings[f"graft_{name}"] = run([
        sys.executable,
        scripts / "postprocess_segment_grafts.py",
        replacement,
        "--output-dir", graft,
        "--short-min-length", 64,
        "--final-min-length", 200,
        "--min-overlap", 31,
        "--overlap-margin", 10,
    ])
    copied = out / f"{name}.fasta"
    shutil.copy2(graft / "primary_contigs.fasta", copied)
    return copied


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
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
        "--bridgeasm", args.bridgeasm,
        "--read1", args.read1,
        "--read2", args.read2,
        "--output", base,
        "--threads", args.threads,
    ])
    current = base / "step6_strain_projection.fasta"

    graph_opt = out / "graph_optimizer"
    timings["graph_optimizer"] = run([
        sys.executable,
        scripts / "graph_path_phaser.py",
        "--gfa", base / "iterative" / "k31_resolve" / "assembly.gfa",
        "-1", args.read1,
        "-2", args.read2,
        "--projection", base / "iterative" / "k21_recall" / "primary_contigs.fasta",
        "--projection", base / "iterative" / "k21_recall" / "haplotigs.fasta",
        "--highk-gfa", base / "iterative" / "k55_resolve" / "assembly.gfa",
        "-o", graph_opt,
        "--anchor-k", 31,
        "--max-context", 6,
        "--dominance", 0.72,
        "--min-direct", 4,
        "--min-length", 200,
    ])

    candidates = {
        "stage1_full_read": graph_opt / "stage1_full_read.fasta",
        "stage2_iterative_projection": graph_opt / "stage2_iterative_projection.fasta",
        "stage3_local_highk": graph_opt / "stage3_local_highk.fasta",
        "stage4_second_pass": graph_opt / "stage4_second_pass.fasta",
    }
    outputs = {"current": str(current)}
    for name, candidate in candidates.items():
        outputs[name] = str(emit(scripts, candidate, current, out, name, timings))

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest = {
        "pipeline": "bridge-na50-segment-graft-v1",
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": usage.ru_maxrss / 1024.0,
        "timings_seconds": timings,
        "outputs": outputs,
    }
    (out / "pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
