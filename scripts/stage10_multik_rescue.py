#!/usr/bin/env python3
"""Stage 10: low-abundance multi-k rescue plus conservative continuity joins.

Stage 8 remains the production backbone. Reads poorly represented by Stage 8
are reassembled independently at k=17/21/25/31. Novel contigs are ranked by
cross-k support so low-abundance sequence can be recovered without simply
relaxing the graph filters. Two rescue candidates are emitted:

* multik_strict: requires cross-k support.
* multik_balanced: also permits very strong single-k contigs under a base cap.

For continuity, a separate pair-backed exact-overlap refinement is emitted.
It calls pair_gap_refine.py with --max-gap 0, so paired reads may choose among
contig-end overlaps but can never create an unknown N-gap scaffold.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import resource
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import low_abundance_rescue as lr


@dataclass
class MultiKCandidate:
    k: int
    source: str
    name: str
    seq: str
    novel_kmers: int
    novel_fraction: float
    cross_k_sources: int
    cross_k_fraction: float
    score: float


def load_raw_candidates(inputs: dict[int, list[Path]], min_length: int) -> list[tuple[int, str, str, str]]:
    raw: list[tuple[int, str, str, str]] = []
    seen: set[tuple[int, str]] = set()
    for k, paths in sorted(inputs.items()):
        for path in paths:
            if not path.exists():
                continue
            for name, seq0 in lr.fasta_records(path):
                seq = lr.canonical(seq0)
                if len(seq) < min_length:
                    continue
                key = (k, seq)
                if key in seen:
                    continue
                seen.add(key)
                raw.append((k, str(path), name, seq))
    return raw


def annotate_multik_candidates(
    raw: list[tuple[int, str, str, str]],
    backbone31: set[str],
    backbone21: set[str],
) -> list[MultiKCandidate]:
    pools: dict[int, set[str]] = defaultdict(set)
    for k, _source, _name, seq in raw:
        pools[k].update(lr.kmers(seq, 21))

    out: list[MultiKCandidate] = []
    for k, source, name, seq in raw:
        all31 = set(lr.kmers(seq, 31))
        if not all31:
            continue
        novel31 = all31 - backbone31
        novel21 = set(lr.kmers(seq, 21)) - backbone21
        supporting_sources = 0
        for other_k, pool in pools.items():
            if other_k == k or not novel21:
                continue
            frac = len(novel21 & pool) / len(novel21)
            if frac >= 0.20:
                supporting_sources += 1
        union_other: set[str] = set()
        for other_k, pool in pools.items():
            if other_k != k:
                union_other.update(pool)
        cross_fraction = (
            len(novel21 & union_other) / len(novel21) if novel21 else 0.0
        )
        novel_fraction = len(novel31) / len(all31)
        consensus_bonus = 1.0 + 0.80 * supporting_sources + 1.25 * cross_fraction
        score = (
            len(novel31)
            * novel_fraction
            * consensus_bonus
            * math.log2(max(2, len(seq)))
        )
        out.append(
            MultiKCandidate(
                k=k,
                source=source,
                name=name,
                seq=seq,
                novel_kmers=len(novel31),
                novel_fraction=novel_fraction,
                cross_k_sources=supporting_sources,
                cross_k_fraction=cross_fraction,
                score=score,
            )
        )
    return out


def select_multik_candidates(
    candidates: list[MultiKCandidate],
    backbone31: set[str],
    *,
    min_novel_kmers: int,
    min_novel_fraction: float,
    min_cross_sources: int,
    min_cross_fraction: float,
    max_total_bases: int,
    max_fraction_per_k: float,
    allow_strong_single_k: bool,
) -> list[MultiKCandidate]:
    eligible: list[MultiKCandidate] = []
    for cand in candidates:
        if cand.novel_kmers < min_novel_kmers or cand.novel_fraction < min_novel_fraction:
            continue
        consensus = (
            cand.cross_k_sources >= min_cross_sources
            and cand.cross_k_fraction >= min_cross_fraction
        )
        strong_single = (
            allow_strong_single_k
            and cand.novel_fraction >= 0.90
            and cand.novel_kmers >= max(100, min_novel_kmers)
            and len(cand.seq) >= 350
        )
        if not consensus and not strong_single:
            continue
        eligible.append(cand)

    eligible.sort(
        key=lambda c: (
            -c.cross_k_sources,
            -c.cross_k_fraction,
            -c.score,
            -c.novel_kmers,
            -len(c.seq),
            c.k,
            c.seq,
        )
    )
    selected: list[MultiKCandidate] = []
    selected_novel: set[str] = set()
    seen_seq: set[str] = set()
    bases_by_k: defaultdict[int, int] = defaultdict(int)
    per_k_cap = max(20_000, int(max_total_bases * max_fraction_per_k))
    total = 0
    for cand in eligible:
        if cand.seq in seen_seq:
            continue
        novel = set(lr.kmers(cand.seq, 31)) - backbone31
        fresh = novel - selected_novel
        if len(fresh) < max(min_novel_kmers, int(0.50 * len(novel))):
            continue
        if bases_by_k[cand.k] + len(cand.seq) > per_k_cap:
            continue
        if total + len(cand.seq) > max_total_bases:
            continue
        selected.append(cand)
        seen_seq.add(cand.seq)
        selected_novel.update(novel)
        bases_by_k[cand.k] += len(cand.seq)
        total += len(cand.seq)
    return selected


def write_metadata(selected: list[MultiKCandidate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "id\tk\tsource\tname\tlength\tnovel_kmers\tnovel_fraction"
            "\tcross_k_sources\tcross_k_fraction\tscore\n"
        )
        for i, c in enumerate(selected, 1):
            handle.write(
                f"{i}\t{c.k}\t{c.source}\t{c.name}\t{len(c.seq)}\t{c.novel_kmers}"
                f"\t{c.novel_fraction:.6f}\t{c.cross_k_sources}"
                f"\t{c.cross_k_fraction:.6f}\t{c.score:.6f}\n"
            )


def run(cmd: list[object], env: dict[str, str] | None = None) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True, env=env)
    return time.monotonic() - started


def assemble_residual_k(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    out: Path,
    k: int,
    threads: int,
) -> float:
    mercy = {17: 72, 21: 56, 25: 44, 31: 32}.get(k, 32)
    env = os.environ.copy()
    if k <= 21:
        env["BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION"] = "0.75"
        env["BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY"] = "35"
        env["BRIDGEASM_MATE_TERMINAL_MERCY_KMERS"] = "96"
    return run(
        [
            bridgeasm,
            "assemble",
            "-1",
            read1,
            "-2",
            read2,
            "-o",
            out,
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
            3,
            "--primary-dominance",
            0.82,
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


def exact_pair_join(
    scripts: Path,
    contigs: Path,
    read1: Path,
    read2: Path,
    output: Path,
    links: Path,
    threads: int,
) -> float:
    return run(
        [
            sys.executable,
            scripts / "pair_gap_refine.py",
            contigs,
            "-1",
            read1,
            "-2",
            read2,
            "-o",
            output,
            "--links",
            links,
            "--threads",
            threads,
            "--min-mapq",
            25,
            "--min-support",
            3,
            "--dominance",
            0.90,
            "--end-window",
            650,
            "--min-overlap",
            31,
            "--max-gap",
            0,
        ]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--ks", default="17,21,25,31")
    ap.add_argument("--pair-stride", type=int, default=3)
    ap.add_argument("--pair-strong-novel", type=float, default=0.15)
    ap.add_argument("--pair-weak-novel", type=float, default=0.45)
    ap.add_argument("--strict-fraction", type=float, default=0.08)
    ap.add_argument("--balanced-fraction", type=float, default=0.12)
    args = ap.parse_args()

    ks = sorted({int(x) for x in args.ks.split(",") if x.strip()})
    if len(ks) < 2 or any(k < 15 or k > 63 for k in ks):
        raise SystemExit("--ks must contain at least two k values in [15,63]")

    started = time.monotonic()
    scripts = Path(__file__).resolve().parent
    out = args.pipeline_dir
    backbone = out / "bridge_backbone.fasta"
    required = [backbone, args.read1, args.read2]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("missing Stage10 inputs: " + ", ".join(missing))

    backbone31, backbone_bases = lr.backbone_kmers(backbone, 31)
    backbone21, _ = lr.backbone_kmers(backbone, 21)
    stage10 = out / "stage10_multik_rescue"
    stage10.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    rare_r1 = stage10 / "rare_R1.fastq.gz"
    rare_r2 = stage10 / "rare_R2.fastq.gz"
    pair_stats = lr.select_residual_pairs(
        args.read1,
        args.read2,
        rare_r1,
        rare_r2,
        backbone31,
        31,
        args.pair_stride,
        args.pair_strong_novel,
        args.pair_weak_novel,
    )

    inputs: dict[int, list[Path]] = {}
    if pair_stats["pairs_kept"] >= 100:
        for k in ks:
            asm = stage10 / f"residual_k{k}"
            timings[f"residual_k{k}_assembly"] = assemble_residual_k(
                args.bridgeasm, rare_r1, rare_r2, asm, k, args.threads
            )
            inputs[k] = [asm / "primary_contigs.fasta", asm / "haplotigs.fasta"]

    raw = load_raw_candidates(inputs, 200)
    candidates = annotate_multik_candidates(raw, backbone31, backbone21)
    strict = select_multik_candidates(
        candidates,
        backbone31,
        min_novel_kmers=64,
        min_novel_fraction=0.70,
        min_cross_sources=1,
        min_cross_fraction=0.30,
        max_total_bases=max(30_000, int(backbone_bases * args.strict_fraction)),
        max_fraction_per_k=0.60,
        allow_strong_single_k=False,
    )
    balanced = select_multik_candidates(
        candidates,
        backbone31,
        min_novel_kmers=40,
        min_novel_fraction=0.50,
        min_cross_sources=1,
        min_cross_fraction=0.20,
        max_total_bases=max(50_000, int(backbone_bases * args.balanced_fraction)),
        max_fraction_per_k=0.60,
        allow_strong_single_k=True,
    )

    strict_add = stage10 / "multik_strict_additions.fasta"
    balanced_add = stage10 / "multik_balanced_additions.fasta"
    lr.write_fasta(
        ((f"multik_strict_{i:06d}_k{c.k}", c.seq) for i, c in enumerate(strict, 1)),
        strict_add,
    )
    lr.write_fasta(
        ((f"multik_balanced_{i:06d}_k{c.k}", c.seq) for i, c in enumerate(balanced, 1)),
        balanced_add,
    )
    write_metadata(strict, stage10 / "multik_strict.tsv")
    write_metadata(balanced, stage10 / "multik_balanced.tsv")

    strict_final = lr.make_union_candidate(
        scripts, backbone, [strict_add], stage10 / "candidate_multik_strict", timings
    )
    balanced_final = lr.make_union_candidate(
        scripts, backbone, [balanced_add], stage10 / "candidate_multik_balanced", timings
    )

    stage8_pair = stage10 / "stage8_pair_exact.fasta"
    strict_pair = stage10 / "multik_strict_pair_exact.fasta"
    balanced_pair = stage10 / "multik_balanced_pair_exact.fasta"
    timings["stage8_pair_exact"] = exact_pair_join(
        scripts,
        backbone,
        args.read1,
        args.read2,
        stage8_pair,
        stage10 / "stage8_pair_exact.links.tsv",
        args.threads,
    )
    timings["strict_pair_exact"] = exact_pair_join(
        scripts,
        strict_final,
        args.read1,
        args.read2,
        strict_pair,
        stage10 / "multik_strict_pair_exact.links.tsv",
        args.threads,
    )
    timings["balanced_pair_exact"] = exact_pair_join(
        scripts,
        balanced_final,
        args.read1,
        args.read2,
        balanced_pair,
        stage10 / "multik_balanced_pair_exact.links.tsv",
        args.threads,
    )

    outputs = {
        "stage8_backbone": str(backbone),
        "stage8_pair_exact": str(stage8_pair),
        "multik_strict": str(strict_final),
        "multik_strict_pair_exact": str(strict_pair),
        "multik_balanced": str(balanced_final),
        "multik_balanced_pair_exact": str(balanced_pair),
    }
    stats = {
        "pipeline": "bridge-stage10-low-abundance-multik-v1",
        "ks": ks,
        "policy": {
            "stage8_replaced": False,
            "cross_k_consensus": True,
            "unknown_gap_scaffolds": False,
            "pair_join": "reciprocal_dominant_pair_backed_exact_overlap",
        },
        "backbone_bases": backbone_bases,
        "residual_pairs": pair_stats,
        "raw_multik_candidates": len(candidates),
        "strict": {
            "selected": len(strict),
            "bases": sum(len(c.seq) for c in strict),
            "by_k": dict(sorted((k, sum(1 for c in strict if c.k == k)) for k in ks)),
        },
        "balanced": {
            "selected": len(balanced),
            "bases": sum(len(c.seq) for c in balanced),
            "by_k": dict(sorted((k, sum(1 for c in balanced if c.k == k)) for k in ks)),
        },
        "timings_seconds": timings,
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0,
        "outputs": outputs,
    }
    (stage10 / "stage10_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    manifest_path = out / "pipeline_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["stage10_low_abundance_multik"] = stats
    manifest.setdefault("outputs", {}).update({f"stage10_{k}": v for k, v in outputs.items()})
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
