#!/usr/bin/env python3
"""Bridge evidence pipeline with physical-fragment-safe multi-k carry-forward.

This is a causal variant of bridge_evidence_pipeline.py.  All downstream stages
are intentionally kept matched.  The only structural change is iterative
handoff:

  previous-k primary/haplotig contigs -> one synthetic physical fragment -> next k

Accepted prior contigs start at 200 bp instead of 500 bp.  All prior sequences
share one physical FASTQ pair and are separated by N^target_k, so any target-k
k-mer receives at most one synthetic fragment support.  It must therefore occur
in at least one real raw fragment to satisfy BridgeAsm's production
min-fragment-support=2 rule.  Unlike sliding virtual pairs, overlapping windows
cannot manufacture multiplicity.

No reference or evaluation metric is used by the assembler.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import sys
import time
from pathlib import Path

import bridge_evidence_pipeline as base


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--prior-min-length", type=int, default=200)
    ap.add_argument("--prior-phred", type=int, default=20)
    args = ap.parse_args()

    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    times: dict[str, float] = {}
    started = time.monotonic()

    step0 = out / "step0_baseline"
    times["step0_baseline"] = base.run(
        [
            sys.executable,
            scripts / "bridge_recovery_pipeline.py",
            "--bridgeasm",
            args.bridgeasm,
            "--read1",
            args.read1,
            "--read2",
            args.read2,
            "--output",
            step0,
            "--threads",
            args.threads,
            "--singleton-fraction",
            0.50,
            "--singleton-quality",
            35,
            "--mate-terminal-mercy",
            96,
            "--stitch-min-overlap",
            31,
            "--stitch-overlap-margin",
            10,
        ]
    )

    iterative = out / "iterative"
    iterative.mkdir(exist_ok=True)
    stages = [
        ("k21_recall", 21, 24),
        ("k31_resolve", 31, 16),
        ("k41_resolve", 41, 12),
        ("k55_resolve", 55, 8),
    ]
    candidates: list[Path] = []
    cur_r1 = args.read1
    cur_r2 = args.read2
    handoff_stats: dict[str, str] = {}

    for idx, (name, k, mercy) in enumerate(stages):
        stage_dir = iterative / name
        env = os.environ.copy()
        if k == 21:
            env["BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION"] = "0.50"
            env["BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY"] = "35"
            env["BRIDGEASM_MATE_TERMINAL_MERCY_KMERS"] = "96"
        else:
            env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION", None)
            env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY", None)
            env.pop("BRIDGEASM_MATE_TERMINAL_MERCY_KMERS", None)

        times[f"step1_{name}"] = base.run(
            base.bridge_cmd(
                args.bridgeasm,
                cur_r1,
                cur_r2,
                stage_dir,
                k,
                mercy,
                args.threads,
            ),
            env=env,
        )
        candidates.append(stage_dir / "primary_contigs.fasta")

        if idx + 1 >= len(stages):
            continue
        next_k = stages[idx + 1][1]
        prior_r1 = iterative / f"{name}.trusted_prior_R1.fastq.gz"
        prior_r2 = iterative / f"{name}.trusted_prior_R2.fastq.gz"
        stats_json = iterative / f"{name}.trusted_prior.json"
        times[f"step1_{name}_trusted_prior"] = base.run(
            [
                sys.executable,
                scripts / "make_trusted_prior_fragment.py",
                stage_dir / "primary_contigs.fasta",
                stage_dir / "haplotigs.fasta",
                "--read1",
                prior_r1,
                "--read2",
                prior_r2,
                "--target-k",
                next_k,
                "--min-length",
                args.prior_min_length,
                "--phred",
                args.prior_phred,
                "--stats-json",
                stats_json,
            ]
        )
        handoff_stats[name] = str(stats_json)
        aug_r1 = iterative / f"{name}.trusted_aug_R1.fastq.gz"
        aug_r2 = iterative / f"{name}.trusted_aug_R2.fastq.gz"
        # Every target-k stage always gets the original physical reads plus one
        # prior fragment.  Synthetic evidence is never recursively copied as
        # ordinary raw reads from an earlier stage.
        base.concat_gzip(args.read1, prior_r1, aug_r1)
        base.concat_gzip(args.read2, prior_r2, aug_r2)
        cur_r1, cur_r2 = aug_r1, aug_r2

    step1_dir = out / "step1_iterative"
    step1 = base.postprocess(scripts, candidates, step1_dir)
    shutil.copy2(step1, out / "step1_iterative.fasta")

    residual = out / "step2_residual_paths.fasta"
    residual_meta = out / "step2_residual_paths.tsv"
    times["step2_residual_extract"] = base.run(
        [
            sys.executable,
            scripts / "residual_path_cover.py",
            iterative / "k31_resolve" / "assembly.gfa",
            iterative / "k41_resolve" / "assembly.gfa",
            iterative / "k55_resolve" / "assembly.gfa",
            "--backbone",
            step1,
            "-o",
            residual,
            "--metadata",
            residual_meta,
            "--secondary-dominance",
            0.35,
            "--extension-dominance",
            0.85,
            "--min-support",
            6,
            "--max-copy",
            2,
            "--novel-k",
            31,
            "--flank",
            120,
            "--max-novel-gap",
            96,
            "--min-novel-kmers",
            4,
            "--min-novel-fraction",
            0.10,
            "--max-patch-length",
            1200,
            "--max-patches",
            300,
            "--max-total-fraction",
            0.05,
        ]
    )
    step2 = base.postprocess(scripts, [step1, residual], out / "step2_residual")
    shutil.copy2(step2, out / "step2_residual.fasta")

    step3 = out / "step3_pairs.fasta"
    times["step3_pair_graph"] = base.run(
        [
            sys.executable,
            scripts / "pair_gap_refine.py",
            step1,
            "-1",
            args.read1,
            "-2",
            args.read2,
            "-o",
            step3,
            "--links",
            out / "step3_pair_links.tsv",
            "--threads",
            args.threads,
            "--min-mapq",
            20,
            "--min-support",
            3,
            "--dominance",
            0.75,
            "--end-window",
            500,
            "--min-overlap",
            31,
            "--max-gap",
            1000,
        ]
    )

    step4 = out / "step4_second_pass.fasta"
    times["step4_second_pass"] = base.run(
        [
            sys.executable,
            scripts / "pair_gap_refine.py",
            step3,
            "-1",
            args.read1,
            "-2",
            args.read2,
            "-o",
            step4,
            "--links",
            out / "step4_pair_links.tsv",
            "--threads",
            args.threads,
            "--min-mapq",
            25,
            "--min-support",
            2,
            "--dominance",
            0.85,
            "--end-window",
            650,
            "--min-overlap",
            31,
            "--max-gap",
            1000,
        ]
    )

    step5 = out / "step5_gapfill.fasta"
    times["step5_gapfill"] = base.run(
        [
            sys.executable,
            scripts / "fill_scaffold_gaps.py",
            step4,
            "-1",
            args.read1,
            "-2",
            args.read2,
            "-o",
            step5,
            "--report",
            out / "step5_gapfill.tsv",
            "--anchor-k",
            31,
            "--local-k",
            21,
            "--flank",
            180,
            "--dominance",
            0.65,
        ]
    )

    projected = out / "step6_projected_strain_paths.fasta"
    projection_map = out / "step6_projection_map.tsv"
    haplotigs = [iterative / name / "haplotigs.fasta" for name, _, _ in stages]
    times["step6_project"] = base.run(
        [
            sys.executable,
            scripts / "strain_projection.py",
            step5,
            *haplotigs,
            "-o",
            projected,
            "--map",
            projection_map,
            "--k",
            31,
            "--projection-k",
            21,
            "--flank-length",
            90,
            "--projection-stride",
            2,
            "--min-flank-hits",
            2,
            "--min-novel-kmers",
            3,
            "--min-novel-fraction",
            0.01,
        ]
    )
    # Step6 is intentionally unchanged.  The experiment isolates the iterative
    # multi-k carry-forward that the breakpoint oracle identified.
    proj_v1 = out / "step6.virtual_R1.fastq.gz"
    proj_v2 = out / "step6.virtual_R2.fastq.gz"
    times["step6_virtualize"] = base.run(
        [
            sys.executable,
            scripts / "make_virtual_pairs.py",
            projected,
            "--read1",
            proj_v1,
            "--read2",
            proj_v2,
            "--read-length",
            101,
            "--insert-size",
            250,
            "--stride",
            120,
            "--min-length",
            250,
        ]
    )
    proj_a1 = out / "step6.aug_R1.fastq.gz"
    proj_a2 = out / "step6.aug_R2.fastq.gz"
    base.concat_gzip(args.read1, proj_v1, proj_a1)
    base.concat_gzip(args.read2, proj_v2, proj_a2)
    projection_asm = out / "step6_projection_k31"
    times["step6_projection_k31"] = base.run(
        base.bridge_cmd(
            args.bridgeasm,
            proj_a1,
            proj_a2,
            projection_asm,
            31,
            16,
            args.threads,
        )
    )
    novel_projection = out / "step6_projection_novel.fasta"
    times["step6_select_novel"] = base.run(
        [
            sys.executable,
            scripts / "strain_projection.py",
            step5,
            projection_asm / "primary_contigs.fasta",
            projection_asm / "haplotigs.fasta",
            "-o",
            novel_projection,
            "--map",
            out / "step6_projection_assembly_map.tsv",
            "--k",
            31,
            "--projection-k",
            21,
            "--flank-length",
            90,
            "--projection-stride",
            2,
            "--min-flank-hits",
            2,
            "--min-novel-kmers",
            3,
            "--min-novel-fraction",
            0.005,
        ]
    )
    step6 = base.postprocess(
        scripts,
        [step5, novel_projection],
        out / "step6_strain_projection",
    )
    shutil.copy2(step6, out / "step6_strain_projection.fasta")

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest = {
        "pipeline": "bridge-evidence-trusted-handoff-v1",
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": usage.ru_maxrss / 1024.0,
        "timings_seconds": times,
        "handoff": {
            "mode": "one_synthetic_physical_fragment_global",
            "min_prior_length": args.prior_min_length,
            "prior_phred": args.prior_phred,
            "target_k_requires_raw_support": True,
            "stats": handoff_stats,
        },
        "promotion_policy": {
            "step1_iterative": "causal_candidate_only",
            "step2_residual": "evaluated_side_candidate_not_promoted",
            "step3_base": "step1_iterative",
        },
        "outputs": {
            "step0_baseline": str(step0 / "primary_contigs.fasta"),
            "step1_iterative": str(out / "step1_iterative.fasta"),
            "step2_residual": str(out / "step2_residual.fasta"),
            "step3_pairs": str(step3),
            "step4_second_pass": str(step4),
            "step5_gapfill": str(step5),
            "step6_strain_projection": str(out / "step6_strain_projection.fasta"),
        },
    }
    (out / "pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
