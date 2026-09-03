#!/usr/bin/env python3
"""Stage 14: amplify the three Stage13 ideas at the scale where they can matter.

A. Read-level superreads: thread Stage10 residual read pairs through the k17
   residual graph, bridge paired anchors only through a unique bounded graph
   path, then extend through unique/context-supported graph edges.  This moves
   superread construction before contig emission instead of overlapping already
   assembled rare contigs.
B. Targeted variable-k: use Stage10 strict rare additions as seed loci, recruit
   full-library read pairs in two taxon-agnostic k-mer rounds, reassemble each
   locus at k17/21/25/31, and permit a higher k to win locally when it reduces
   branching or raises N50 while retaining enough sequence.
C. Long-component flow: enumerate 0.3--5 kb source->sink tangles on the k31
   graph, fit several paths simultaneously by non-negative coverage
   decomposition, and retain physically/read-thread supported alternative paths.

All emitted sequence still comes from an existing graph or a local assembly.
No reference is used.  Each method is tested independently and through a
conservative combined candidate on top of Stage10 strict.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import resource
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import adaptive_k_local_v2 as ak
import graph_path_phaser as gp
import low_abundance_rescue as lr
import stage10_multik_rescue as s10
import stage13_three_methods as s13


@dataclass
class Tangle:
    source: str
    sink: str
    paths: list[list[str]]
    min_bp: int
    max_bp: int


@dataclass
class LocalAssembly:
    seed_id: int
    k: int
    output_dir: Path
    bases: int
    n50: int
    contigs: int
    branches: int


def run(cmd: list[object], env: dict[str, str] | None = None) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True, env=env)
    return time.monotonic() - started


def write_fasta(records: Iterable[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name} len={len(seq)}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")


def path_bp(graph: gp.Graph, path: list[str]) -> int:
    if not path:
        return 0
    return len(gp.path_sequence(path, graph))


def path_key(graph: gp.Graph, path: list[str]) -> tuple[str, ...]:
    forward = tuple(path)
    reverse = tuple(graph.rev.get(uid, uid) for uid in reversed(path))
    return min(forward, reverse)


def edge_score(graph: gp.Graph, left: str, right: str) -> float:
    ev = graph.edge.get((left, right), gp.EdgeEvidence())
    return ev.direct * 2.0 + ev.gapped + ev.pairs * 3.0


def unique_bounded_bridge(
    graph: gp.Graph,
    source: str,
    target: str,
    *,
    max_edges: int = 14,
    max_bp: int = 900,
    max_states: int = 2500,
) -> list[str] | None:
    """Return a unique source->target graph path under edge/bp bounds."""
    if source == target:
        return [source]
    stack: list[tuple[list[str], int]] = [([source], len(graph.seqs[source]))]
    solutions: list[list[str]] = []
    states = 0
    while stack and states < max_states:
        path, bases = stack.pop()
        states += 1
        if len(path) - 1 >= max_edges:
            continue
        node = path[-1]
        for child in reversed(graph.out.get(node, [])):
            if child in path or graph.rev.get(child, child) in path:
                continue
            added = max(1, len(graph.seqs[child]) - graph.k)
            new_bases = bases + added
            if new_bases > max_bp:
                continue
            path2 = path + [child]
            if child == target:
                solutions.append(path2)
                if len(solutions) > 1:
                    return None
                continue
            stack.append((path2, new_bases))
    return solutions[0] if len(solutions) == 1 else None


def _choose_context_extension(
    graph: gp.Graph,
    history: list[str],
    raw_ctx: Counter[tuple[str, ...]],
    *,
    forward: bool,
) -> str | None:
    current = history[-1] if forward else history[0]
    adjacent = graph.out.get(current, []) if forward else graph.inc.get(current, [])
    used = set(history) | {graph.rev.get(uid, uid) for uid in history}
    candidates = [uid for uid in adjacent if uid not in used and graph.rev.get(uid, uid) not in used]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    scored: list[tuple[float, int, int, str]] = []
    for uid in candidates:
        support, span = gp.context_strength(raw_ctx, history, uid, forward)
        left, right = (current, uid) if forward else (uid, current)
        score = support * max(2, span) * 8.0 + edge_score(graph, left, right)
        scored.append((score, support, span, uid))
    scored.sort(key=lambda x: (-x[0], -x[2], -x[1], x[3]))
    best = scored[0]
    second = scored[1][0] if len(scored) > 1 else 0.0
    if best[1] >= 2 and best[2] >= 3 and best[0] >= max(1.0, second * 1.35):
        return best[3]
    current_uid = current
    left, right = (current_uid, best[3]) if forward else (best[3], current_uid)
    ev = graph.edge.get((left, right), gp.EdgeEvidence())
    if (ev.pairs >= 2 or ev.direct >= 4) and best[0] >= max(4.0, second * 1.60):
        return best[3]
    return None


def extend_threaded_path(
    graph: gp.Graph,
    path0: list[str],
    raw_ctx: Counter[tuple[str, ...]],
    *,
    max_nodes: int = 48,
    max_bp: int = 3000,
) -> list[str]:
    path = list(path0)
    if not path:
        return path
    for forward in (False, True):
        while len(path) < max_nodes and path_bp(graph, path) < max_bp:
            nxt = _choose_context_extension(graph, path, raw_ctx, forward=forward)
            if nxt is None:
                break
            if forward:
                path.append(nxt)
            else:
                path.insert(0, nxt)
    return path


def longest_thread(seq: str, graph: gp.Graph, index: gp.KmerIndex) -> list[str] | None:
    segments = gp.thread_sequence(seq, graph, index, None, max_bridge_edges=3)
    if not segments:
        return None
    return max(segments, key=lambda path: (path_bp(graph, path), len(path), tuple(path)))


def collect_read_superreads(
    graph: gp.Graph,
    index: gp.KmerIndex,
    read1: Path,
    read2: Path,
    raw_ctx: Counter[tuple[str, ...]],
) -> tuple[list[str], dict[str, int]]:
    seen_paths: set[tuple[str, ...]] = set()
    sequences: dict[str, str] = {}
    pair_closures = 0
    threaded_pairs = 0
    unique_extensions = 0
    bridge_cache: dict[tuple[str, str], list[str] | None] = {}

    left_iter = gp.read_fastq(read1)
    right_iter = gp.read_fastq(read2)
    for left, right in zip(left_iter, right_iter):
        r1 = left[1]
        r2 = right[1]
        pair_had_thread = False
        orientations = ((r1, gp.rc(r2)), (r2, gp.rc(r1)))
        for left_seq, right_seq in orientations:
            lpath = longest_thread(left_seq, graph, index)
            rpath = longest_thread(right_seq, graph, index)
            if lpath is None and rpath is None:
                continue
            pair_had_thread = True
            for path0 in (lpath, rpath):
                if path0 is None:
                    continue
                extended = extend_threaded_path(graph, path0, raw_ctx)
                if len(extended) > len(path0):
                    unique_extensions += 1
                key = path_key(graph, extended)
                if key in seen_paths:
                    continue
                seen_paths.add(key)
                seq = lr.canonical(gp.path_sequence(extended, graph))
                if len(seq) >= 220:
                    sequences[seq] = seq
            if lpath is None or rpath is None:
                continue
            cache_key = (lpath[-1], rpath[0])
            if cache_key not in bridge_cache:
                bridge_cache[cache_key] = unique_bounded_bridge(
                    graph, lpath[-1], rpath[0], max_edges=14, max_bp=900
                )
            bridge = bridge_cache[cache_key]
            if bridge is None:
                continue
            combined = list(lpath)
            combined.extend(bridge[1:-1])
            combined.extend(rpath)
            collapsed: list[str] = []
            for uid in combined:
                if not collapsed or uid != collapsed[-1]:
                    collapsed.append(uid)
            combined = extend_threaded_path(graph, collapsed, raw_ctx)
            key = path_key(graph, combined)
            if key in seen_paths:
                continue
            seen_paths.add(key)
            seq = lr.canonical(gp.path_sequence(combined, graph))
            if len(seq) >= 260:
                sequences[seq] = seq
                pair_closures += 1
        if pair_had_thread:
            threaded_pairs += 1
    ordered = sorted(sequences, key=lambda seq: (-len(seq), seq))
    return ordered, {
        "threaded_pairs": threaded_pairs,
        "pair_closures": pair_closures,
        "unique_extensions": unique_extensions,
        "raw_superreads": len(ordered),
        "raw_superread_bases": sum(map(len, ordered)),
        "bridge_cache_entries": len(bridge_cache),
    }


def select_bridge_sequences(
    seqs: Iterable[tuple[str, str]],
    strict_baseline: Path,
    pools: dict[int, set[str]],
    *,
    min_length: int,
    min_novel_kmers: int,
    min_novel_fraction: float,
    min_cross_sources: int,
    min_cross_fraction: float,
    max_total_bases: int,
) -> list[s13.SeqEvidence]:
    base31, _ = lr.backbone_kmers(strict_baseline, 31)
    base21, _ = lr.backbone_kmers(strict_baseline, 21)
    items: list[s13.SeqEvidence] = []
    for name, seq in seqs:
        item = s13.sequence_evidence(name, seq, base31, base21, pools)
        if item is not None:
            items.append(item)
    return s13.select_sequence_evidence(
        items,
        base31,
        min_length=min_length,
        min_novel_kmers=min_novel_kmers,
        min_novel_fraction=min_novel_fraction,
        min_cross_sources=min_cross_sources,
        min_cross_fraction=min_cross_fraction,
        min_fresh_fraction=0.45,
        max_total_bases=max_total_bases,
    )


def make_bridge_candidate(
    scripts: Path,
    baseline: Path,
    additions: list[Path],
    outdir: Path,
    timings: dict[str, float],
    *,
    min_overlap: int = 81,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    union = outdir / "union.fasta"
    noncontained = outdir / "noncontained.fasta"
    final = outdir / "primary_contigs.fasta"
    timings[f"merge_{outdir.name}"] = run(
        [sys.executable, scripts / "merge_fasta_unique.py", union, baseline, *additions, "--min-length", 200]
    )
    timings[f"contain_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "filter_contained_fasta.py",
            union,
            noncontained,
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
        ]
    )
    timings[f"stitch_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "stitch_exact_overlaps.py",
            final,
            noncontained,
            "--min-overlap",
            min_overlap,
            "--overlap-margin",
            20,
            "--seed-length",
            31,
            "--max-seed-occurrences",
            16,
            "--min-length",
            200,
        ]
    )
    return final


def build_read_superread_candidate(
    scripts: Path,
    stage10: Path,
    strict_baseline: Path,
    pools: dict[int, set[str]],
    outdir: Path,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    graph = gp.Graph.from_gfa(stage10 / "residual_k17" / "assembly.gfa")
    index = gp.KmerIndex(graph, 17)
    rare1 = stage10 / "rare_R1.fastq.gz"
    rare2 = stage10 / "rare_R2.fastq.gz"
    started = time.monotonic()
    raw_ctx, ctx_stats = gp.collect_read_contexts(graph, index, rare1, rare2, None, 12)
    seqs, thread_stats = collect_read_superreads(graph, index, rare1, rare2, raw_ctx)
    timings["read_superread_threading"] = time.monotonic() - started
    _, baseline_bases = lr.backbone_kmers(strict_baseline, 31)
    selected = select_bridge_sequences(
        ((f"threaded_{i:07d}", seq) for i, seq in enumerate(seqs, 1)),
        strict_baseline,
        pools,
        min_length=260,
        min_novel_kmers=24,
        min_novel_fraction=0.20,
        min_cross_sources=1,
        min_cross_fraction=0.30,
        max_total_bases=max(80_000, int(baseline_bases * 0.12)),
    )
    additions = outdir / "read_superread_additions.fasta"
    write_fasta(((f"stage14_superread_{i:06d}", item.seq) for i, item in enumerate(selected, 1)), additions)
    s13.write_evidence(selected, outdir / "read_superread_selected.tsv")
    final = make_bridge_candidate(
        scripts, strict_baseline, [additions], outdir / "candidate_read_superread", timings, min_overlap=81
    )
    return final, additions, {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "contexts": ctx_stats,
        **thread_stats,
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
    }


def seed_signature_index(
    seed_fasta: Path,
    backbone: Path,
    *,
    k: int = 21,
    max_seeds: int = 36,
    min_signature_kmers: int = 12,
) -> tuple[list[tuple[str, str]], dict[str, int], dict[int, set[str]]]:
    backbone_kmers, _ = lr.backbone_kmers(backbone, k)
    seeds = [(name, lr.canonical(seq)) for name, seq in lr.fasta_records(seed_fasta)]
    memberships: dict[str, set[int]] = defaultdict(set)
    per_seed: dict[int, set[str]] = {}
    for sid, (_name, seq) in enumerate(seeds):
        novel = set(lr.kmers(seq, k)) - backbone_kmers
        per_seed[sid] = novel
        for mer in novel:
            memberships[mer].add(sid)
    unique_by_seed: dict[int, set[str]] = defaultdict(set)
    for mer, ids in memberships.items():
        if len(ids) == 1:
            unique_by_seed[next(iter(ids))].add(mer)
    ranked = sorted(
        (sid for sid in range(len(seeds)) if len(unique_by_seed[sid]) >= min_signature_kmers),
        key=lambda sid: (-len(unique_by_seed[sid]), -len(seeds[sid][1]), seeds[sid][1]),
    )[:max_seeds]
    remap = {old: new for new, old in enumerate(ranked)}
    chosen = [seeds[old] for old in ranked]
    signatures: dict[str, int] = {}
    chosen_sets: dict[int, set[str]] = {}
    for old, new in remap.items():
        chosen_sets[new] = unique_by_seed[old]
        for mer in unique_by_seed[old]:
            signatures[mer] = new
    return chosen, signatures, chosen_sets


def pair_scores(seq1: str, seq2: str, signatures: dict[str, int], k: int, stride: int) -> Counter[int]:
    scores: Counter[int] = Counter()
    for seq in (seq1, seq2):
        for mer in lr.kmers(seq, k, stride):
            sid = signatures.get(mer)
            if sid is not None:
                scores[sid] += 1
    return scores


def assign_seed(scores: Counter[int], *, min_hits: int, margin: int) -> int | None:
    ranked = scores.most_common(2)
    if not ranked or ranked[0][1] < min_hits:
        return None
    if len(ranked) > 1 and ranked[0][1] - ranked[1][1] < margin:
        return None
    return ranked[0][0]


def recruit_seed_pairs(
    read1: Path,
    read2: Path,
    seeds: list[tuple[str, str]],
    initial_signatures: dict[str, int],
    backbone: Path,
    *,
    initial_k: int = 21,
    expand_k: int = 17,
) -> tuple[dict[int, set[int]], dict[str, int]]:
    assignments: dict[int, set[int]] = defaultdict(set)
    initial_pairs = 0
    for idx, (left, right) in enumerate(zip(ak.fastq_records(read1), ak.fastq_records(read2))):
        scores = pair_scores(left[1], right[1], initial_signatures, initial_k, 2)
        sid = assign_seed(scores, min_hits=2, margin=1)
        if sid is not None:
            assignments[sid].add(idx)
            initial_pairs += 1

    backbone17, _ = lr.backbone_kmers(backbone, expand_k)
    memberships: dict[str, set[int]] = defaultdict(set)
    for idx, (left, right) in enumerate(zip(ak.fastq_records(read1), ak.fastq_records(read2))):
        hit_seeds = [sid for sid, ids in assignments.items() if idx in ids]
        if not hit_seeds:
            continue
        sid = hit_seeds[0]
        local: set[str] = set()
        local.update(lr.kmers(left[1], expand_k, 2))
        local.update(lr.kmers(right[1], expand_k, 2))
        for mer in local:
            if mer not in backbone17:
                memberships[mer].add(sid)
    expanded = {
        mer: next(iter(ids))
        for mer, ids in memberships.items()
        if len(ids) == 1
    }

    second_pairs = 0
    for idx, (left, right) in enumerate(zip(ak.fastq_records(read1), ak.fastq_records(read2))):
        if any(idx in ids for ids in assignments.values()):
            continue
        scores = pair_scores(left[1], right[1], expanded, expand_k, 2)
        sid = assign_seed(scores, min_hits=3, margin=2)
        if sid is not None:
            assignments[sid].add(idx)
            second_pairs += 1
    return assignments, {
        "seed_loci": len(seeds),
        "initial_signatures": len(initial_signatures),
        "expanded_signatures": len(expanded),
        "initial_pairs": initial_pairs,
        "second_round_pairs": second_pairs,
        "total_assigned_pairs": sum(len(ids) for ids in assignments.values()),
    }


def write_seed_fastqs(
    read1: Path,
    read2: Path,
    assignments: dict[int, set[int]],
    outdir: Path,
) -> dict[int, tuple[Path, Path]]:
    outdir.mkdir(parents=True, exist_ok=True)
    handles: dict[int, tuple[object, object]] = {}
    paths: dict[int, tuple[Path, Path]] = {}
    for sid in assignments:
        d = outdir / f"seed_{sid:03d}" / "reads"
        d.mkdir(parents=True, exist_ok=True)
        p1, p2 = d / "R1.fastq.gz", d / "R2.fastq.gz"
        handles[sid] = (gzip.open(p1, "wt"), gzip.open(p2, "wt"))
        paths[sid] = (p1, p2)
    for idx, (left, right) in enumerate(zip(ak.fastq_records(read1), ak.fastq_records(read2))):
        for sid, ids in assignments.items():
            if idx not in ids:
                continue
            out1, out2 = handles[sid]
            ak.write_record(out1, left)
            ak.write_record(out2, right)
    for out1, out2 in handles.values():
        out1.close()
        out2.close()
    return paths


def assemble_local_seed(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    outdir: Path,
    seed_id: int,
    k: int,
    threads: int,
) -> LocalAssembly:
    s10.assemble_residual_k(bridgeasm, read1, read2, outdir, k, threads)
    profile = json.loads((outdir / "run_profile.json").read_text())
    return LocalAssembly(
        seed_id=seed_id,
        k=k,
        output_dir=outdir,
        bases=int(profile.get("primary_bases", 0)),
        n50=int(profile.get("primary_n50", 0)),
        contigs=int(profile.get("primary_contigs", 0)),
        branches=ak.gfa_branch_count(outdir / "assembly.gfa"),
    )


def choose_variable_k(candidates: list[LocalAssembly], base_k: int = 17) -> tuple[LocalAssembly | None, str]:
    valid = [c for c in candidates if c.bases > 0]
    if not valid:
        return None, "no_output"
    base = next((c for c in candidates if c.k == base_k), None)
    if base is None or base.bases == 0:
        best = max(valid, key=lambda c: (c.n50, -c.branches, c.bases, c.k))
        return best, "base_empty_choose_best"
    promoted: list[LocalAssembly] = []
    for c in valid:
        if c.k <= base_k:
            continue
        retain = c.bases / max(1, base.bases)
        branch_ratio = c.branches / max(1, base.branches)
        n50_ratio = c.n50 / max(1, base.n50)
        contig_ratio = c.contigs / max(1, base.contigs)
        if retain < 0.50:
            continue
        if branch_ratio <= 0.75 or n50_ratio >= 1.15 or (contig_ratio <= 0.75 and n50_ratio >= 1.0):
            promoted.append(c)
    if not promoted:
        return base, "keep_k17"
    def score(c: LocalAssembly) -> tuple[float, int, int, int]:
        retain = c.bases / max(1, base.bases)
        branch_gain = 1.0 - c.branches / max(1, base.branches)
        n50_gain = math.log2(max(1.0, c.n50 / max(1, base.n50)))
        return (1.8 * branch_gain + 1.5 * n50_gain + 0.4 * min(1.2, retain), c.n50, -c.branches, c.k)
    return max(promoted, key=score), "promote_higher_k"


def build_variable_k_candidate(
    scripts: Path,
    bridgeasm: Path,
    stage10: Path,
    strict_baseline: Path,
    backbone: Path,
    pools: dict[int, set[str]],
    read1: Path,
    read2: Path,
    outdir: Path,
    threads: int,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    seed_fasta = stage10 / "multik_strict_additions.fasta"
    seeds, signatures, _ = seed_signature_index(seed_fasta, backbone)
    recruitment_started = time.monotonic()
    assignments, recruit_stats = recruit_seed_pairs(read1, read2, seeds, signatures, backbone)
    seed_fastqs = write_seed_fastqs(read1, read2, assignments, outdir / "targeted_reads")
    timings["variable_k_recruitment"] = time.monotonic() - recruitment_started

    choices: list[LocalAssembly] = []
    all_records: list[tuple[str, str]] = []
    promoted = 0
    assembled_loci = 0
    choice_rows: list[tuple[object, ...]] = []
    for sid, (r1, r2) in sorted(seed_fastqs.items()):
        pair_count = len(assignments[sid])
        if pair_count < 24:
            choice_rows.append((sid, pair_count, 0, 0, 0, 0, False, "too_few_pairs"))
            continue
        assembled_loci += 1
        local_candidates: list[LocalAssembly] = []
        for k in (17, 21, 25, 31):
            asm = outdir / "targeted_reads" / f"seed_{sid:03d}" / f"k{k}"
            local_candidates.append(assemble_local_seed(bridgeasm, r1, r2, asm, sid, k, threads))
        chosen, reason = choose_variable_k(local_candidates)
        if chosen is None:
            choice_rows.append((sid, pair_count, 0, 0, 0, 0, False, reason))
            continue
        if chosen.k > 17:
            promoted += 1
        choices.append(chosen)
        choice_rows.append((sid, pair_count, chosen.k, chosen.bases, chosen.n50, chosen.branches, chosen.k > 17, reason))
        for fasta in (chosen.output_dir / "primary_contigs.fasta", chosen.output_dir / "haplotigs.fasta"):
            if not fasta.exists():
                continue
            for name, seq in lr.fasta_records(fasta):
                all_records.append((f"seed{sid}_k{chosen.k}_{name}", seq))
    with (outdir / "variable_k_choices.tsv").open("w") as handle:
        handle.write("seed\tpairs\tchosen_k\tbases\tn50\tbranches\tpromoted\treason\n")
        for row in choice_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    _, baseline_bases = lr.backbone_kmers(strict_baseline, 31)
    selected = select_bridge_sequences(
        all_records,
        strict_baseline,
        pools,
        min_length=240,
        min_novel_kmers=20,
        min_novel_fraction=0.18,
        min_cross_sources=1,
        min_cross_fraction=0.28,
        max_total_bases=max(100_000, int(baseline_bases * 0.15)),
    )
    additions = outdir / "variable_k_additions.fasta"
    write_fasta(((f"stage14_variable_k_{i:06d}", item.seq) for i, item in enumerate(selected, 1)), additions)
    s13.write_evidence(selected, outdir / "variable_k_selected.tsv")
    final = make_bridge_candidate(
        scripts, strict_baseline, [additions], outdir / "candidate_variable_k", timings, min_overlap=81
    )
    return final, additions, {
        **recruit_stats,
        "assembled_loci": assembled_loci,
        "higher_k_promotions": promoted,
        "selected_loci": len(choices),
        "candidate_records": len(all_records),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
    }


def _enumerate_tangle_paths(
    graph: gp.Graph,
    source: str,
    *,
    max_depth: int,
    max_bp: int,
    max_branch: int,
    max_states: int,
) -> dict[str, list[list[str]]]:
    by_sink: dict[str, list[list[str]]] = defaultdict(list)
    states = 0
    for child in graph.out.get(source, []):
        stack: list[tuple[list[str], int]] = [([source, child], path_bp(graph, [source, child]))]
        while stack and states < max_states:
            path, bases = stack.pop()
            states += 1
            node = path[-1]
            if len(path) >= 3 and len(graph.inc.get(node, [])) > 1:
                by_sink[node].append(path)
            if len(path) - 1 >= max_depth or bases >= max_bp:
                continue
            children = graph.out.get(node, [])
            if not children or len(children) > max_branch:
                continue
            for nxt in reversed(children):
                if nxt in path or graph.rev.get(nxt, nxt) in path:
                    continue
                added = max(1, len(graph.seqs[nxt]) - graph.k)
                if bases + added > max_bp:
                    continue
                stack.append((path + [nxt], bases + added))
    return by_sink


def discover_long_tangles(
    graph: gp.Graph,
    *,
    min_bp: int = 300,
    max_bp: int = 5000,
    max_depth: int = 40,
    max_paths: int = 8,
    max_branch: int = 4,
    max_tangles: int = 220,
) -> list[Tangle]:
    tangles: list[Tangle] = []
    seen: set[tuple[str, str]] = set()
    def source_priority(uid: str) -> tuple[float, float, int, str]:
        covs = sorted(
            (max(0.0, graph.coverage.get(child, 0.0)) for child in graph.out.get(uid, [])),
            reverse=True,
        )
        total = sum(covs)
        minor_share = (covs[1] / total) if len(covs) > 1 and total > 0 else 0.0
        mixture_score = min(minor_share, 0.5) * max(0.0, 0.55 - minor_share)
        return (-mixture_score, -graph.coverage.get(uid, 0.0), -len(graph.seqs[uid]), uid)

    sources = sorted(
        (uid for uid in graph.seqs if 2 <= len(graph.out.get(uid, [])) <= max_branch),
        key=source_priority,
    )[:1800]
    for source in sources:
        sinks = _enumerate_tangle_paths(
            graph,
            source,
            max_depth=max_depth,
            max_bp=max_bp,
            max_branch=max_branch,
            max_states=12000,
        )
        ranked: list[tuple[float, str, list[list[str]], int, int]] = []
        for sink, paths0 in sinks.items():
            unique = {tuple(path): path for path in paths0}
            paths = list(unique.values())
            if not (2 <= len(paths) <= max_paths):
                continue
            if len({path[1] for path in paths if len(path) > 1}) < 2:
                continue
            lengths = [path_bp(graph, path) for path in paths]
            if min(lengths) < min_bp or max(lengths) > max_bp:
                continue
            score = min(lengths) + 0.25 * sum(lengths) / len(lengths)
            ranked.append((score, sink, paths, min(lengths), max(lengths)))
        if not ranked:
            continue
        _score, sink, paths, lo, hi = max(ranked, key=lambda x: (x[0], -len(x[2]), x[1]))
        key = min((source, sink), (graph.rev.get(sink, sink), graph.rev.get(source, source)))
        if key in seen:
            continue
        seen.add(key)
        paths.sort(key=lambda p: (path_bp(graph, p), tuple(p)))
        tangles.append(Tangle(source, sink, paths, lo, hi))
        if len(tangles) >= max_tangles:
            break
    return tangles


def path_thread_fraction(counter: Counter[tuple[str, ...]], path: list[str]) -> tuple[int, float, int]:
    edges = max(0, len(path) - 1)
    supported = sum(counter.get((a, b), 0) > 0 for a, b in zip(path, path[1:]))
    max_support = max((counter.get((a, b), 0) for a, b in zip(path, path[1:])), default=0)
    return supported, supported / max(1, edges), max_support


def flow_sequence_evidence(
    name: str,
    seq0: str,
    baseline31: set[str],
    pools: dict[int, set[str]],
) -> tuple[s13.SeqEvidence | None, int, float]:
    seq = lr.canonical(seq0)
    all31 = set(lr.kmers(seq, 31))
    if not all31:
        return None, 0, 0.0
    novel31 = all31 - baseline31
    novel_fraction = len(novel31) / len(all31)
    dummy21: set[str] = set()
    item = s13.sequence_evidence(name, seq, baseline31, dummy21, pools)
    return item, len(novel31), novel_fraction


def build_long_flow_candidate(
    scripts: Path,
    strict_baseline: Path,
    target_gfa: Path,
    read1: Path,
    read2: Path,
    pools: dict[int, set[str]],
    outdir: Path,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    graph = gp.Graph.from_gfa(target_gfa)
    index = gp.KmerIndex(graph, 31)
    started = time.monotonic()
    raw_ctx, raw_stats = gp.collect_read_contexts(graph, index, read1, read2, None, 12)
    tangles = discover_long_tangles(graph)
    timings["long_flow_threading_discovery"] = time.monotonic() - started
    baseline31, baseline_bases = lr.backbone_kmers(strict_baseline, 31)
    lengths = {uid: len(seq) for uid, seq in graph.seqs.items()}
    eligible: list[s13.SeqEvidence] = []
    rows: list[tuple[object, ...]] = []
    raw_paths = 0
    for tid, tangle in enumerate(tangles, 1):
        abundances, rmse = s13.projected_nnls(tangle.paths, graph.coverage, lengths, iterations=100, l1=0.03)
        total = sum(abundances)
        for pid, (path, abundance) in enumerate(zip(tangle.paths, abundances), 1):
            raw_paths += 1
            share = abundance / total if total > 0 else 0.0
            physical, edge_count, pair_support = s13.path_physical_evidence(graph, path)
            physical_fraction = physical / max(1, edge_count)
            thread_edges, thread_fraction, max_thread = path_thread_fraction(raw_ctx, path)
            seq = gp.path_sequence(path, graph)
            item, novel_kmers, novel_fraction = flow_sequence_evidence(
                f"tangle{tid}_path{pid}", seq, baseline31, pools
            )
            cross_sources = item.cross_k_sources if item is not None else 0
            cross_fraction = item.cross_k_fraction if item is not None else 0.0
            rows.append(
                (
                    tid,
                    tangle.source,
                    tangle.sink,
                    pid,
                    len(path),
                    len(seq),
                    abundance,
                    share,
                    rmse,
                    physical,
                    edge_count,
                    physical_fraction,
                    pair_support,
                    thread_edges,
                    thread_fraction,
                    max_thread,
                    novel_kmers,
                    novel_fraction,
                    cross_sources,
                    cross_fraction,
                )
            )
            if len(seq) < 300 or novel_kmers < 24 or novel_fraction < 0.05:
                continue
            if abundance < 0.50 or share < 0.05:
                continue
            if physical_fraction < 0.80 or thread_fraction < 0.30:
                continue
            if max_thread < 1 and pair_support < 2:
                continue
            if item is None:
                score = novel_kmers * novel_fraction * math.log2(max(2, len(seq)))
                item = s13.SeqEvidence(
                    name=f"tangle{tid}_path{pid}",
                    seq=lr.canonical(seq),
                    novel_kmers=novel_kmers,
                    novel_fraction=novel_fraction,
                    cross_k_sources=0,
                    cross_k_fraction=0.0,
                    score=score,
                )
            eligible.append(item)
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "long_flow_paths.tsv").open("w") as handle:
        handle.write(
            "tangle\tsource\tsink\tpath_index\tnodes\tlength\tabundance\tshare\trmse"
            "\tphysical_edges\tedge_count\tphysical_fraction\tpair_support\tthread_edges"
            "\tthread_fraction\tmax_thread\tnovel_kmers\tnovel_fraction\tcross_k_sources"
            "\tcross_k_fraction\n"
        )
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")

    eligible.sort(
        key=lambda item: (-item.cross_k_sources, -item.cross_k_fraction, -item.novel_kmers, -len(item.seq), item.seq)
    )
    selected: list[s13.SeqEvidence] = []
    selected_novel: set[str] = set()
    total_bases = 0
    max_bases = max(120_000, int(baseline_bases * 0.18))
    for item in eligible:
        novel = set(lr.kmers(item.seq, 31)) - baseline31
        fresh = novel - selected_novel
        if len(fresh) < max(20, math.ceil(0.35 * len(novel))):
            continue
        if total_bases + len(item.seq) > max_bases:
            continue
        selected.append(item)
        selected_novel.update(novel)
        total_bases += len(item.seq)
    additions = outdir / "long_flow_additions.fasta"
    write_fasta(((f"stage14_long_flow_{i:06d}", item.seq) for i, item in enumerate(selected, 1)), additions)
    s13.write_evidence(selected, outdir / "long_flow_selected.tsv")
    final = make_bridge_candidate(
        scripts, strict_baseline, [additions], outdir / "candidate_long_flow", timings, min_overlap=81
    )
    return final, additions, {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "tangles": len(tangles),
        "raw_flow_paths": raw_paths,
        "eligible_paths": len(eligible),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "read_threading": raw_stats,
    }


def build_combined_candidate(
    scripts: Path,
    strict_baseline: Path,
    additions: list[Path],
    outdir: Path,
    timings: dict[str, float],
) -> tuple[Path, dict[str, int]]:
    final = make_bridge_candidate(
        scripts, strict_baseline, additions, outdir / "candidate_combined", timings, min_overlap=81
    )
    records = 0
    bases = 0
    for path in additions:
        if not path.exists():
            continue
        for _name, seq in lr.fasta_records(path):
            records += 1
            bases += len(seq)
    return final, {"input_records": records, "input_bases": bases}


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
    inputs = s13.residual_inputs(stage10, ks)
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
        stage10 / "multik_strict_additions.fasta",
    ]
    for paths in inputs.values():
        required.extend(paths)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage14 inputs: " + ", ".join(missing))

    stage14 = out / "stage14_amplified_methods"
    stage14.mkdir(parents=True, exist_ok=True)
    pools = s13.cross_k_pools(inputs)
    timings: dict[str, float] = {}

    super_final, super_add, super_stats = build_read_superread_candidate(
        scripts, stage10, strict_baseline, pools, stage14 / "read_superread", timings
    )
    var_final, var_add, var_stats = build_variable_k_candidate(
        scripts,
        args.bridgeasm,
        stage10,
        strict_baseline,
        backbone,
        pools,
        args.read1,
        args.read2,
        stage14 / "variable_k",
        args.threads,
        timings,
    )
    flow_final, flow_add, flow_stats = build_long_flow_candidate(
        scripts,
        strict_baseline,
        target_gfa,
        args.read1,
        args.read2,
        pools,
        stage14 / "long_flow",
        timings,
    )
    combined_final, combined_stats = build_combined_candidate(
        scripts,
        strict_baseline,
        [super_add, var_add, flow_add],
        stage14 / "combined",
        timings,
    )

    stats = {
        "pipeline": "bridge-stage14-amplified-methods-v1",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "superread_source": "Stage10 residual read pairs threaded on residual k17 graph",
            "variable_k_seed": "Stage10 strict rare additions; two-round full-library recruitment",
            "flow_scale": "k31 source-to-sink tangles 300-5000 bp",
            "sequence_join": "reciprocal-best exact overlap >=81 bp",
        },
        "methods": {
            "read_superread": super_stats,
            "variable_k": var_stats,
            "long_flow": flow_stats,
            "combined": combined_stats,
        },
        "outputs": {
            "read_superread": str(super_final),
            "variable_k": str(var_final),
            "long_flow": str(flow_final),
            "combined": str(combined_final),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (stage14 / "stage14_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
