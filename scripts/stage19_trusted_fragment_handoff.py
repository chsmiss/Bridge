#!/usr/bin/env python3
"""Stage19: carry validated rare priors as exactly one synthetic fragment.

This is a reference-free causal experiment for the production multi-k handoff.
All trusted contigs at one stage are concatenated into one synthetic R1 record,
separated by N^k.  R2 is a single N.  BridgeAsm counts fragment support once per
physical pair after deduplicating k-mers, so the entire prior contributes at
most +1 fragment support to any target-k k-mer, globally.  A k-mer therefore
still needs occurrence in a real raw fragment to pass the production
min-fragment-support=2 gate.

Unlike sliding virtual pairs this cannot self-promote a prior k-mer by overlap,
and unlike a 500-bp cutoff it can carry validated 200--499 bp rare paths.  Q20
is used for the synthetic record: the prior supplies structural support but
cannot raise a low-quality raw observation above the production mean-Q20 gate.

After each k31/k41/k55 assembly only records connected to a locus-unique Stage10
seed signature are carried forward.  Final additions also require recurrence at
another target k through the Stage16 cross-k evidence routine.
"""
from __future__ import annotations

import argparse
import gzip
import json
import resource
import time
from collections import Counter
from pathlib import Path

import low_abundance_rescue as lr
import stage14_amplified_methods as s14
import stage16_root_cause as s16
import stage18_short_handoff as s18


def write_trusted_fragment(
    source_fasta: Path,
    out1: Path,
    out2: Path,
    *,
    target_k: int,
    synthetic_phred: int = 20,
) -> dict[str, int | float]:
    if not (0 <= synthetic_phred <= 40):
        raise ValueError("synthetic_phred must be in 0..=40")
    sequences: list[str] = []
    seen: set[str] = set()
    source_bases = 0
    for _name, seq0 in lr.fasta_records(source_fasta):
        seq = seq0.upper()
        if len(seq) < target_k or "N" in seq:
            continue
        canonical = lr.canonical(seq)
        if canonical in seen:
            continue
        seen.add(canonical)
        sequences.append(seq)
        source_bases += len(seq)

    out1.parent.mkdir(parents=True, exist_ok=True)
    separator = "N" * target_k
    joined = separator.join(sequences)
    qchar = chr(33 + synthetic_phred)
    with gzip.open(out1, "wt", compresslevel=3) as left, gzip.open(
        out2, "wt", compresslevel=3
    ) as right:
        if joined:
            left.write(f"@trusted_prior_k{target_k}/1\n{joined}\n+\n{qchar * len(joined)}\n")
            right.write(f"@trusted_prior_k{target_k}/2\nN\n+\n!\n")

    distinct_target_kmers: set[str] = set()
    target_kmer_observations = 0
    for seq in sequences:
        kmers = list(lr.kmers(seq, target_k))
        target_kmer_observations += len(kmers)
        distinct_target_kmers.update(kmers)
    return {
        "source_records": len(sequences),
        "source_bases": source_bases,
        "synthetic_fragments": int(bool(joined)),
        "synthetic_r1_bases": len(joined),
        "separator_bases": max(0, len(sequences) - 1) * target_k,
        "target_kmer_observations": target_kmer_observations,
        "distinct_target_kmers": len(distinct_target_kmers),
        "max_prior_fragment_support": int(bool(joined)),
        "synthetic_phred": synthetic_phred,
    }


def build_trusted_fragment_handoff(
    scripts: Path,
    bridgeasm: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    backbone: Path,
    read1: Path,
    read2: Path,
    threads: int,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    stage10 = pipeline_dir / "stage10_multik_rescue"
    seed_fasta = stage10 / "multik_strict_additions.fasta"
    root = pipeline_dir / "stage19_trusted_fragment_handoff"
    root.mkdir(parents=True, exist_ok=True)

    seeds, signatures, _sets, signature_stats = s16.build_all_seed_signatures(
        seed_fasta, backbone, k=21, min_signature_kmers=4
    )
    source = seed_fasta
    assembly_inputs: dict[int, list[Path]] = {}
    stage_stats: dict[str, object] = {}

    for k in (31, 41, 55):
        stage = root / f"k{k}"
        stage.mkdir(parents=True, exist_ok=True)
        prior1 = stage / "trusted_prior_R1.fastq.gz"
        prior2 = stage / "trusted_prior_R2.fastq.gz"
        prior_stats = write_trusted_fragment(
            source,
            prior1,
            prior2,
            target_k=k,
            synthetic_phred=20,
        )
        aug1 = stage / "aug_R1.fastq.gz"
        aug2 = stage / "aug_R2.fastq.gz"
        s18.concat_gzip([read1, prior1], aug1)
        s18.concat_gzip([read2, prior2], aug2)
        asm = stage / "assembly"
        timings[f"trusted_fragment_k{k}"] = s18.assemble_target_k(
            bridgeasm, aug1, aug2, asm, k, threads
        )
        connected_path = stage / "seed_connected.fasta"
        connected, connected_stats = s18.filter_seed_connected(
            [asm / "primary_contigs.fasta", asm / "haplotigs.fasta"],
            signatures,
            connected_path,
        )
        assembly_inputs[k] = [connected_path]
        stage_stats[f"k{k}"] = {
            "prior": prior_stats,
            "assembly_profile": s18.profile(asm / "run_profile.json"),
            "connected": connected_stats,
            "connected_lengths": {
                "lt250": sum(len(item.seq) < 250 for item in connected),
                "250_to_499": sum(250 <= len(item.seq) < 500 for item in connected),
                "ge500": sum(len(item.seq) >= 500 for item in connected),
                "max": max((len(item.seq) for item in connected), default=0),
            },
        }
        source = connected_path

    evidence = s16.local_evidence(
        assembly_inputs,
        seeds,
        strict_baseline,
        Counter(),
    )
    selected = s16.select_local_evidence(evidence, strict_baseline)
    s16.write_local_evidence(evidence, root / "trusted_fragment_evidence.tsv")
    s16.write_local_evidence(selected, root / "trusted_fragment_selected.tsv")
    additions = root / "trusted_fragment_additions.fasta"
    s14.write_fasta(
        (
            (f"stage19_seed{item.seed_id}_k{item.k}_{idx:06d}", item.seq)
            for idx, item in enumerate(selected, 1)
        ),
        additions,
    )
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        root / "candidate_trusted_fragment",
        timings,
        min_overlap=81,
    )
    return final, additions, {
        **signature_stats,
        "stages": stage_stats,
        "evidence_records": len(evidence),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "selected_seed_loci": len({item.seed_id for item in selected}),
        "selected_fresh31": sum(item.fresh31 for item in selected),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    scripts = Path(__file__).resolve().parent
    pipeline = args.pipeline_dir
    stage10 = pipeline / "stage10_multik_rescue"
    strict_baseline = stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    backbone = pipeline / "bridge_backbone.fasta"
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
        raise SystemExit("missing Stage19 inputs: " + ", ".join(missing))

    started = time.monotonic()
    timings: dict[str, float] = {}
    final, additions, method_stats = build_trusted_fragment_handoff(
        scripts,
        args.bridgeasm,
        pipeline,
        strict_baseline,
        backbone,
        args.read1,
        args.read2,
        args.threads,
        timings,
    )
    root = pipeline / "stage19_trusted_fragment_handoff"
    stats = {
        "pipeline": "bridge-stage19-trusted-fragment-handoff-v1",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "metric_targets": False,
            "prior_source": "Stage10 strict cross-k validated rare additions",
            "handoff": "all priors in one synthetic physical fragment separated by N^k",
            "max_prior_fragment_support_per_target_kmer": 1,
            "synthetic_phred": 20,
            "raw_requirement": "production min-fragment-support=2 therefore at least one real fragment is still required",
            "next_stage_filter": "locus-unique Stage10 seed signatures",
            "final_validation": "seed-connected + Stage10-novel + cross-target-k recurrence",
        },
        "method": method_stats,
        "outputs": {"final": str(final), "additions": str(additions)},
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (root / "stage19_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
