#!/usr/bin/env python3
"""NA50-focused cumulative graph optimization pipeline.

Starts from the current evidence-six-step assembly, then replaces the dominant
k31 path extraction with four cumulative graph-only refinements:
  1 full-read path phasing
  2 low-k graph path projection
  3 high-k evidence only at ambiguous target-graph junctions
  4 context-aware second-pass read threading and graph re-resolution

Each optimized graph assembly is merged with the current recovery output only
through containment/exact-overlap post-processing so low-abundance sequence is
not discarded while the dominant backbone can become more contiguous.
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


def postprocess(scripts: Path, inputs: list[Path], outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    merged = outdir / "union.fasta"
    filtered = outdir / "noncontained.fasta"
    final = outdir / "primary_contigs.fasta"
    run([sys.executable, scripts / "merge_fasta_unique.py", merged, *inputs, "--min-length", 200])
    run([
        sys.executable,
        scripts / "filter_contained_fasta.py",
        merged,
        filtered,
        "--min-length", 200,
        "--seed-k", 21,
        "--window", 12,
        "--candidate-minimizers", 16,
        "--removed-tsv", outdir / "contained_removed.tsv",
        "--stats-json", outdir / "containment_stats.json",
    ])
    run([
        sys.executable,
        scripts / "stitch_exact_overlaps.py",
        final,
        filtered,
        "--min-overlap", 31,
        "--overlap-margin", 10,
        "--seed-length", 31,
        "--max-seed-occurrences", 64,
        "--min-length", 200,
    ])
    return final


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

    stage_files = {
        "stage1_full_read": graph_opt / "stage1_full_read.fasta",
        "stage2_iterative_projection": graph_opt / "stage2_iterative_projection.fasta",
        "stage3_local_highk": graph_opt / "stage3_local_highk.fasta",
        "stage4_second_pass": graph_opt / "stage4_second_pass.fasta",
    }
    outputs: dict[str, str] = {"current": str(current)}
    for name, candidate in stage_files.items():
        stage_dir = out / name
        t0 = time.monotonic()
        final = postprocess(scripts, [current, candidate], stage_dir)
        timings[f"postprocess_{name}"] = time.monotonic() - t0
        copied = out / f"{name}.fasta"
        shutil.copy2(final, copied)
        outputs[name] = str(copied)

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest = {
        "pipeline": "bridge-na50-graph-v1",
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": usage.ru_maxrss / 1024.0,
        "timings_seconds": timings,
        "outputs": outputs,
    }
    (out / "pipeline_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
