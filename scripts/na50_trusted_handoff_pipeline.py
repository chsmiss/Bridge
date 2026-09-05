#!/usr/bin/env python3
"""Matched Stage8-style graph optimization for trusted multi-k handoff.

The graph/phasing/repeat/postprocess configuration is copied from
na50_repeat_pipeline.py.  Only the underlying current_pipeline producer is
swapped to bridge_evidence_trusted_handoff_pipeline.py, keeping the causal
comparison focused on multi-k carry-forward.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import na50_repeat_pipeline as base


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--segment-anchor-bases", type=int, default=31)
    ap.add_argument("--prior-min-length", type=int, default=200)
    ap.add_argument("--prior-phred", type=int, default=20)
    args = ap.parse_args()
    if args.segment_anchor_bases < 31:
        raise SystemExit("segment-anchor-bases must be >=31")

    scripts = Path(__file__).resolve().parent
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}
    started = time.monotonic()

    current_pipeline = out / "current_pipeline"
    timings["current_pipeline"] = base.run(
        [
            sys.executable,
            scripts / "bridge_evidence_trusted_handoff_pipeline.py",
            "--bridgeasm",
            args.bridgeasm,
            "--read1",
            args.read1,
            "--read2",
            args.read2,
            "--output",
            current_pipeline,
            "--threads",
            args.threads,
            "--prior-min-length",
            args.prior_min_length,
            "--prior-phred",
            args.prior_phred,
        ]
    )
    current = current_pipeline / "step6_strain_projection.fasta"
    target_gfa = current_pipeline / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = current_pipeline / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = current_pipeline / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = current_pipeline / "iterative" / "k55_resolve" / "assembly.gfa"

    graph_opt = out / "graph_optimizer"
    timings["graph_optimizer"] = base.run(
        [
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
        ]
    )

    repeat_opt = out / "repeat_optimizer"
    timings["repeat_optimizer"] = base.run(
        [
            sys.executable,
            scripts / "repeat_graph_optimizer_v2.py",
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
        ]
    )

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
            base.emit_stage(
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
        "pipeline": "bridge-na50-trusted-handoff-v1",
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": usage.ru_maxrss / 1024.0,
        "segment_anchor_bases": args.segment_anchor_bases,
        "prior_min_length": args.prior_min_length,
        "prior_phred": args.prior_phred,
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
