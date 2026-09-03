#!/usr/bin/env python3
"""Stage 13: three cross-domain assembly prototypes on top of Stage10.

The production baseline stays Stage8 + Stage10 strict low-abundance rescue.
Stage13 evaluates three orthogonal ideas before any core Rust rewrite:

A. Super-read closure: relaxed-but-cross-k-supported rare contigs are joined only
   through reciprocal-best exact suffix/prefix overlaps. Only genuinely longer
   merged records with strong cross-k support are promoted.
B. Local multiplex DBG: bounded branch neighborhoods from the residual k17 graph
   are routed their own read pairs and locally reassembled at k17/21/25/31.
   Output still has to pass Stage10-style novelty and cross-k gates.
C. Local multi-path flow: small k31 bubbles are decomposed into multiple paths at
   once. Non-negative path abundances are fitted to node coverage; low-abundance
   paths are rescued only when read-thread/edge evidence and multi-k evidence
   agree. Coverage is therefore a component-level constraint, not a greedy
   junction tie breaker.

Each method emits an independent candidate plus a conservative combined
candidate so MetaQUAST can reveal which abstraction actually moves GF/NA50.
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
import stage10_multik_rescue as s10


@dataclass
class SeqEvidence:
    name: str
    seq: str
    novel_kmers: int
    novel_fraction: float
    cross_k_sources: int
    cross_k_fraction: float
    score: float


@dataclass
class Bubble:
    source: str
    sink: str
    paths: list[list[str]]


@dataclass
class FlowCandidate:
    bubble: int
    path_index: int
    path: list[str]
    abundance: float
    abundance_share: float
    thread_support: int
    physical_edges: int
    edge_count: int
    pair_support: int
    evidence: SeqEvidence


def run(cmd: list[object]) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True)
    return time.monotonic() - started


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name} len={len(seq)}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")


def residual_inputs(stage10: Path, ks: list[int]) -> dict[int, list[Path]]:
    inputs: dict[int, list[Path]] = {}
    for k in ks:
        directory = stage10 / f"residual_k{k}"
        inputs[k] = [directory / "primary_contigs.fasta", directory / "haplotigs.fasta"]
    return inputs


def cross_k_pools(inputs: dict[int, list[Path]], kmer: int = 21) -> dict[int, set[str]]:
    pools: dict[int, set[str]] = defaultdict(set)
    for k, paths in inputs.items():
        for path in paths:
            if not path.exists():
                continue
            for _name, seq in lr.fasta_records(path):
                pools[k].update(lr.kmers(seq, kmer))
    return dict(pools)


def sequence_evidence(
    name: str,
    seq0: str,
    baseline31: set[str],
    baseline21: set[str],
    pools: dict[int, set[str]],
) -> SeqEvidence | None:
    seq = lr.canonical(seq0)
    all31 = set(lr.kmers(seq, 31))
    if not all31:
        return None
    novel31 = all31 - baseline31
    novel21 = set(lr.kmers(seq, 21)) - baseline21
    if not novel21:
        return None
    source_fractions = [
        len(novel21 & pool) / len(novel21)
        for pool in pools.values()
        if pool
    ]
    cross_sources = sum(frac >= 0.20 for frac in source_fractions)
    union: set[str] = set()
    for pool in pools.values():
        union.update(pool)
    cross_fraction = len(novel21 & union) / len(novel21)
    novel_fraction = len(novel31) / len(all31)
    score = (
        len(novel31)
        * novel_fraction
        * (1.0 + cross_sources + 1.5 * cross_fraction)
        * math.log2(max(2, len(seq)))
    )
    return SeqEvidence(
        name=name,
        seq=seq,
        novel_kmers=len(novel31),
        novel_fraction=novel_fraction,
        cross_k_sources=cross_sources,
        cross_k_fraction=cross_fraction,
        score=score,
    )


def select_sequence_evidence(
    candidates: list[SeqEvidence],
    baseline31: set[str],
    *,
    min_length: int,
    min_novel_kmers: int,
    min_novel_fraction: float,
    min_cross_sources: int,
    min_cross_fraction: float,
    min_fresh_fraction: float,
    max_total_bases: int,
) -> list[SeqEvidence]:
    eligible = [
        item
        for item in candidates
        if len(item.seq) >= min_length
        and item.novel_kmers >= min_novel_kmers
        and item.novel_fraction >= min_novel_fraction
        and item.cross_k_sources >= min_cross_sources
        and item.cross_k_fraction >= min_cross_fraction
    ]
    eligible.sort(
        key=lambda item: (
            -item.cross_k_sources,
            -item.cross_k_fraction,
            -item.novel_fraction,
            -item.score,
            -len(item.seq),
            item.seq,
        )
    )
    selected: list[SeqEvidence] = []
    selected_novel: set[str] = set()
    seen: set[str] = set()
    total = 0
    for item in eligible:
        if item.seq in seen:
            continue
        novel = set(lr.kmers(item.seq, 31)) - baseline31
        fresh = novel - selected_novel
        required = max(min_novel_kmers, math.ceil(len(novel) * min_fresh_fraction))
        if len(fresh) < required:
            continue
        if total + len(item.seq) > max_total_bases:
            continue
        selected.append(item)
        selected_novel.update(novel)
        seen.add(item.seq)
        total += len(item.seq)
    return selected


def write_evidence(items: list[SeqEvidence], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "id\tname\tlength\tnovel_kmers\tnovel_fraction\tcross_k_sources"
            "\tcross_k_fraction\tscore\n"
        )
        for i, item in enumerate(items, 1):
            handle.write(
                f"{i}\t{item.name}\t{len(item.seq)}\t{item.novel_kmers}\t"
                f"{item.novel_fraction:.6f}\t{item.cross_k_sources}\t"
                f"{item.cross_k_fraction:.6f}\t{item.score:.6f}\n"
            )


def build_superread_candidate(
    scripts: Path,
    stage10: Path,
    strict_baseline: Path,
    backbone: Path,
    inputs: dict[int, list[Path]],
    pools: dict[int, set[str]],
    outdir: Path,
    timings: dict[str, float],
) -> tuple[Path, dict[str, object], Path]:
    backbone31, backbone_bases = lr.backbone_kmers(backbone, 31)
    backbone21, _ = lr.backbone_kmers(backbone, 21)
    raw = s10.load_raw_candidates(inputs, 200)
    annotated = s10.annotate_multik_candidates(raw, backbone31, backbone21)
    sources = s10.select_multik_candidates(
        annotated,
        backbone31,
        min_novel_kmers=32,
        min_novel_fraction=0.55,
        min_cross_sources=1,
        min_cross_fraction=0.25,
        max_total_bases=max(80_000, int(backbone_bases * 0.15)),
        max_fraction_per_k=0.70,
        allow_strong_single_k=False,
    )
    source_fasta = outdir / "superread_sources.fasta"
    write_fasta(
        [(f"source_{i:06d}_k{item.k}", item.seq) for i, item in enumerate(sources, 1)],
        source_fasta,
    )
    stitched = outdir / "superreads.stitched.fasta"
    timings["superread_exact_closure"] = run(
        [
            sys.executable,
            scripts / "stitch_exact_overlaps.py",
            stitched,
            source_fasta,
            "--min-overlap",
            45,
            "--overlap-margin",
            8,
            "--seed-length",
            21,
            "--min-length",
            200,
        ]
    )
    source_seqs = {lr.canonical(item.seq) for item in sources}
    evidence: list[SeqEvidence] = []
    merged_records = 0
    for name, seq in lr.fasta_records(stitched):
        can = lr.canonical(seq)
        if can in source_seqs:
            continue
        merged_records += 1
        item = sequence_evidence(name, can, backbone31, backbone21, pools)
        if item is not None:
            evidence.append(item)
    selected = select_sequence_evidence(
        evidence,
        backbone31,
        min_length=300,
        min_novel_kmers=48,
        min_novel_fraction=0.65,
        min_cross_sources=2,
        min_cross_fraction=0.55,
        min_fresh_fraction=0.55,
        max_total_bases=max(50_000, int(backbone_bases * 0.10)),
    )
    additions = outdir / "superread_additions.fasta"
    write_fasta(
        [(f"stage13_superread_{i:06d}", item.seq) for i, item in enumerate(selected, 1)],
        additions,
    )
    write_evidence(selected, outdir / "superread_selected.tsv")
    final = lr.make_union_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_superread",
        timings,
    )
    stats = {
        "source_records": len(sources),
        "source_bases": sum(len(item.seq) for item in sources),
        "merged_records": merged_records,
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
    }
    return final, stats, additions


def build_multiplex_candidate(
    scripts: Path,
    bridgeasm: Path,
    stage10: Path,
    strict_baseline: Path,
    backbone: Path,
    pools: dict[int, set[str]],
    outdir: Path,
    threads: int,
    timings: dict[str, float],
) -> tuple[Path, dict[str, object], Path]:
    base_dir = stage10 / "residual_k17"
    rare1 = stage10 / "rare_R1.fastq.gz"
    rare2 = stage10 / "rare_R2.fastq.gz"
    local = outdir / "multiplex_local"
    timings["multiplex_local_reassembly"] = run(
        [
            sys.executable,
            scripts / "adaptive_k_local_v2.py",
            "--bridgeasm",
            bridgeasm,
            "--read1",
            rare1,
            "--read2",
            rare2,
            "--base-dir",
            base_dir,
            "--output",
            local,
            "--base-k",
            17,
            "--candidate-k",
            "21,25,31",
            "--seed-k",
            15,
            "--radius",
            2,
            "--max-nodes",
            48,
            "--max-neighborhoods",
            24,
            "--min-branch-nodes",
            2,
            "--min-pairs",
            20,
            "--min-candidate-base-fraction",
            0.85,
            "--branch-fraction",
            0.85,
            "--n50-gain",
            1.10,
            "--threads",
            threads,
        ]
    )
    backbone31, backbone_bases = lr.backbone_kmers(backbone, 31)
    backbone21, _ = lr.backbone_kmers(backbone, 21)
    evidence: list[SeqEvidence] = []
    adaptive = local / "adaptive_contigs.fasta"
    if adaptive.exists():
        for name, seq in lr.fasta_records(adaptive):
            item = sequence_evidence(name, seq, backbone31, backbone21, pools)
            if item is not None:
                evidence.append(item)
    selected = select_sequence_evidence(
        evidence,
        backbone31,
        min_length=250,
        min_novel_kmers=48,
        min_novel_fraction=0.70,
        min_cross_sources=2,
        min_cross_fraction=0.45,
        min_fresh_fraction=0.55,
        max_total_bases=max(50_000, int(backbone_bases * 0.10)),
    )
    additions = outdir / "multiplex_additions.fasta"
    write_fasta(
        [(f"stage13_multiplex_{i:06d}", item.seq) for i, item in enumerate(selected, 1)],
        additions,
    )
    write_evidence(selected, outdir / "multiplex_selected.tsv")
    final = lr.make_union_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_multiplex",
        timings,
    )
    summary = {}
    summary_path = local / "summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
    stats = {
        "local_summary": summary,
        "candidate_records": len(evidence),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
    }
    return final, stats, additions


def _enumerate_branch_paths(
    graph: gp.Graph,
    source: str,
    max_depth: int,
    max_branch: int,
    max_states: int,
) -> dict[str, list[list[str]]]:
    by_sink: dict[str, list[list[str]]] = defaultdict(list)
    states = 0
    for child in graph.out.get(source, []):
        stack: list[list[str]] = [[source, child]]
        while stack and states < max_states:
            path = stack.pop()
            states += 1
            node = path[-1]
            depth = len(path) - 1
            if depth >= 2 and len(graph.inc.get(node, [])) > 1:
                by_sink[node].append(path)
            if depth >= max_depth:
                continue
            children = graph.out.get(node, [])
            if not children or len(children) > max_branch:
                continue
            for nxt in reversed(children):
                if nxt in path or graph.rev.get(nxt, nxt) in path:
                    continue
                stack.append(path + [nxt])
    return by_sink


def discover_bubbles(
    graph: gp.Graph,
    *,
    max_depth: int,
    max_paths: int,
    max_branch: int,
    max_bubbles: int,
) -> list[Bubble]:
    bubbles: list[Bubble] = []
    seen_keys: set[tuple[str, str]] = set()
    sources = sorted(
        (
            uid
            for uid in graph.seqs
            if 2 <= len(graph.out.get(uid, [])) <= max_branch
        ),
        key=lambda uid: (
            -len(graph.out.get(uid, [])),
            -graph.coverage.get(uid, 0.0),
            uid,
        ),
    )
    for source in sources:
        sinks = _enumerate_branch_paths(
            graph, source, max_depth, max_branch, max_states=4000
        )
        ranked: list[tuple[tuple[int, int, str], str, list[list[str]]]] = []
        for sink, paths0 in sinks.items():
            unique: dict[tuple[str, ...], list[str]] = {tuple(path): path for path in paths0}
            paths = list(unique.values())
            first_branches = {path[1] for path in paths if len(path) >= 2}
            if len(first_branches) < 2 or len(paths) < 2 or len(paths) > max_paths:
                continue
            max_len = max(len(path) for path in paths)
            if max_len > max_depth + 1:
                continue
            ranked.append(((max_len, sum(map(len, paths)), sink), sink, paths))
        if not ranked:
            continue
        _rank, sink, paths = min(ranked, key=lambda item: item[0])
        rev_source = graph.rev.get(source, source)
        rev_sink = graph.rev.get(sink, sink)
        key = min((source, sink), (rev_sink, rev_source))
        if key in seen_keys:
            continue
        seen_keys.add(key)
        paths.sort(key=lambda path: (len(path), tuple(path)))
        bubbles.append(Bubble(source=source, sink=sink, paths=paths))
        if len(bubbles) >= max_bubbles:
            break
    return bubbles


def projected_nnls(
    paths: list[list[str]],
    coverage: dict[str, float],
    lengths: dict[str, int],
    *,
    iterations: int = 80,
    l1: float = 0.05,
) -> tuple[list[float], float]:
    nodes = sorted({uid for path in paths for uid in path[1:-1]})
    if not nodes:
        return [0.0] * len(paths), 0.0
    memberships = [set(path[1:-1]) for path in paths]
    abundances = []
    for members in memberships:
        vals = [max(0.0, coverage.get(uid, 0.0)) for uid in members]
        abundances.append(min(vals) if vals else 0.0)
    for _ in range(iterations):
        max_delta = 0.0
        for j, members in enumerate(memberships):
            if not members:
                continue
            numerator = 0.0
            denominator = 0.0
            for uid in members:
                weight = max(1.0, float(lengths.get(uid, 1)))
                other = 0.0
                for k, other_members in enumerate(memberships):
                    if k != j and uid in other_members:
                        other += abundances[k]
                numerator += weight * (max(0.0, coverage.get(uid, 0.0)) - other)
                denominator += weight
            new_value = max(0.0, (numerator - l1) / max(denominator, 1.0))
            max_delta = max(max_delta, abs(new_value - abundances[j]))
            abundances[j] = new_value
        if max_delta < 1e-4:
            break
    squared = 0.0
    weight_sum = 0.0
    for uid in nodes:
        weight = max(1.0, float(lengths.get(uid, 1)))
        fitted = sum(
            abundances[j]
            for j, members in enumerate(memberships)
            if uid in members
        )
        residual = max(0.0, coverage.get(uid, 0.0)) - fitted
        squared += weight * residual * residual
        weight_sum += weight
    rmse = math.sqrt(squared / max(weight_sum, 1.0))
    return abundances, rmse


def path_thread_support(counter: Counter[tuple[str, ...]], path: list[str]) -> int:
    best = 0
    for size in range(2, min(6, len(path)) + 1):
        for start in range(0, len(path) - size + 1):
            best = max(best, counter.get(tuple(path[start : start + size]), 0))
    return best


def path_physical_evidence(graph: gp.Graph, path: list[str]) -> tuple[int, int, int]:
    edge_count = max(0, len(path) - 1)
    physical_edges = 0
    pair_support = 0
    for left, right in zip(path, path[1:]):
        ev = graph.edge.get((left, right), gp.EdgeEvidence())
        support = ev.direct + ev.gapped + ev.pairs
        physical_edges += int(support > 0)
        pair_support += ev.pairs
    return physical_edges, edge_count, pair_support


def build_flow_candidate(
    stage10: Path,
    strict_baseline: Path,
    backbone: Path,
    target_gfa: Path,
    read1: Path,
    read2: Path,
    pools: dict[int, set[str]],
    outdir: Path,
    scripts: Path,
    timings: dict[str, float],
) -> tuple[Path, dict[str, object], Path]:
    graph = gp.Graph.from_gfa(target_gfa)
    index = gp.KmerIndex(graph, 31)
    context_started = time.monotonic()
    raw_ctx, raw_stats = gp.collect_read_contexts(graph, index, read1, read2, None, 10)
    timings["flow_read_threading"] = time.monotonic() - context_started
    bubbles = discover_bubbles(
        graph,
        max_depth=8,
        max_paths=8,
        max_branch=4,
        max_bubbles=600,
    )
    baseline31, baseline_bases = lr.backbone_kmers(strict_baseline, 31)
    baseline21, _ = lr.backbone_kmers(strict_baseline, 21)
    raw_candidates: list[FlowCandidate] = []
    report_rows: list[tuple[object, ...]] = []
    lengths = {uid: len(seq) for uid, seq in graph.seqs.items()}
    for bid, bubble in enumerate(bubbles, 1):
        abundances, rmse = projected_nnls(bubble.paths, graph.coverage, lengths)
        total_abundance = sum(abundances)
        for pid, (path, abundance) in enumerate(zip(bubble.paths, abundances), 1):
            share = abundance / total_abundance if total_abundance > 0 else 0.0
            thread_support = path_thread_support(raw_ctx, path)
            physical_edges, edge_count, pair_support = path_physical_evidence(graph, path)
            internal = path[1:-1]
            if not internal:
                continue
            seq = gp.path_sequence(internal, graph)
            item = sequence_evidence(
                f"bubble{bid}_path{pid}", seq, baseline31, baseline21, pools
            )
            if item is None:
                continue
            candidate = FlowCandidate(
                bubble=bid,
                path_index=pid,
                path=path,
                abundance=abundance,
                abundance_share=share,
                thread_support=thread_support,
                physical_edges=physical_edges,
                edge_count=edge_count,
                pair_support=pair_support,
                evidence=item,
            )
            raw_candidates.append(candidate)
            report_rows.append(
                (
                    bid,
                    bubble.source,
                    bubble.sink,
                    pid,
                    ",".join(path),
                    abundance,
                    share,
                    rmse,
                    thread_support,
                    physical_edges,
                    edge_count,
                    pair_support,
                    len(seq),
                    item.novel_kmers,
                    item.novel_fraction,
                    item.cross_k_sources,
                    item.cross_k_fraction,
                )
            )
    with (outdir / "flow_paths.tsv").open("w") as handle:
        handle.write(
            "bubble\tsource\tsink\tpath_index\tpath\tabundance\tshare\trmse"
            "\tthread_support\tphysical_edges\tedge_count\tpair_support\tlength"
            "\tnovel_kmers\tnovel_fraction\tcross_k_sources\tcross_k_fraction\n"
        )
        for row in report_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    eligible: list[SeqEvidence] = []
    for cand in raw_candidates:
        edge_fraction = cand.physical_edges / max(1, cand.edge_count)
        evidence_ok = (
            cand.thread_support >= 1
            or (cand.pair_support >= 2 and edge_fraction >= 0.75)
        )
        if not evidence_ok:
            continue
        if cand.abundance < 0.75 or cand.abundance_share < 0.08:
            continue
        item = cand.evidence
        if item.cross_k_sources < 1 or item.cross_k_fraction < 0.35:
            continue
        eligible.append(item)
    selected = select_sequence_evidence(
        eligible,
        baseline31,
        min_length=200,
        min_novel_kmers=32,
        min_novel_fraction=0.40,
        min_cross_sources=1,
        min_cross_fraction=0.35,
        min_fresh_fraction=0.55,
        max_total_bases=max(40_000, int(baseline_bases * 0.08)),
    )
    additions = outdir / "flow_additions.fasta"
    write_fasta(
        [(f"stage13_flow_{i:06d}", item.seq) for i, item in enumerate(selected, 1)],
        additions,
    )
    write_evidence(selected, outdir / "flow_selected.tsv")
    final = lr.make_union_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_flow",
        timings,
    )
    stats = {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "bubbles": len(bubbles),
        "raw_flow_paths": len(raw_candidates),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "read_threading": raw_stats,
    }
    return final, stats, additions


def build_combined_candidate(
    scripts: Path,
    strict_baseline: Path,
    pools: dict[int, set[str]],
    addition_files: list[Path],
    outdir: Path,
    timings: dict[str, float],
) -> tuple[Path, dict[str, object]]:
    baseline31, baseline_bases = lr.backbone_kmers(strict_baseline, 31)
    baseline21, _ = lr.backbone_kmers(strict_baseline, 21)
    evidence: list[SeqEvidence] = []
    for path in addition_files:
        if not path.exists():
            continue
        for name, seq in lr.fasta_records(path):
            item = sequence_evidence(name, seq, baseline31, baseline21, pools)
            if item is not None:
                evidence.append(item)
    selected = select_sequence_evidence(
        evidence,
        baseline31,
        min_length=200,
        min_novel_kmers=32,
        min_novel_fraction=0.40,
        min_cross_sources=1,
        min_cross_fraction=0.35,
        min_fresh_fraction=0.60,
        max_total_bases=max(70_000, int(baseline_bases * 0.12)),
    )
    additions = outdir / "combined_additions.fasta"
    write_fasta(
        [(f"stage13_combined_{i:06d}", item.seq) for i, item in enumerate(selected, 1)],
        additions,
    )
    write_evidence(selected, outdir / "combined_selected.tsv")
    final = lr.make_union_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_combined",
        timings,
    )
    return final, {
        "input_records": len(evidence),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
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
    backbone = out / "bridge_backbone.fasta"
    strict_baseline = stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    target_gfa = out / "current_pipeline" / "iterative" / "k31_resolve" / "assembly.gfa"
    ks = sorted({int(value) for value in args.ks.split(",") if value.strip()})
    inputs = residual_inputs(stage10, ks)
    required = [
        args.bridgeasm,
        args.read1,
        args.read2,
        backbone,
        strict_baseline,
        target_gfa,
        stage10 / "rare_R1.fastq.gz",
        stage10 / "rare_R2.fastq.gz",
        stage10 / "residual_k17" / "assembly.gfa",
    ]
    for paths in inputs.values():
        required.extend(paths)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage13 inputs: " + ", ".join(missing))

    stage13 = out / "stage13_three_methods"
    stage13.mkdir(parents=True, exist_ok=True)
    pools = cross_k_pools(inputs)
    timings: dict[str, float] = {}

    superread_final, superread_stats, superread_add = build_superread_candidate(
        scripts,
        stage10,
        strict_baseline,
        backbone,
        inputs,
        pools,
        stage13 / "superread",
        timings,
    )
    multiplex_final, multiplex_stats, multiplex_add = build_multiplex_candidate(
        scripts,
        args.bridgeasm,
        stage10,
        strict_baseline,
        backbone,
        pools,
        stage13 / "multiplex",
        args.threads,
        timings,
    )
    flow_final, flow_stats, flow_add = build_flow_candidate(
        stage10,
        strict_baseline,
        backbone,
        target_gfa,
        args.read1,
        args.read2,
        pools,
        stage13 / "flow",
        scripts,
        timings,
    )
    combined_final, combined_stats = build_combined_candidate(
        scripts,
        strict_baseline,
        pools,
        [superread_add, multiplex_add, flow_add],
        stage13 / "combined",
        timings,
    )

    stats = {
        "pipeline": "bridge-stage13-three-methods-v1",
        "baseline": str(strict_baseline),
        "methods": {
            "superread": superread_stats,
            "multiplex": multiplex_stats,
            "flow": flow_stats,
            "combined": combined_stats,
        },
        "outputs": {
            "superread": str(superread_final),
            "multiplex": str(multiplex_final),
            "flow": str(flow_final),
            "combined": str(combined_final),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (stage13 / "stage13_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
