#!/usr/bin/env python3
"""Stage 17: exact-preserving graph threading plus iterative locus frontier recovery.

Stage16 showed two concrete failure modes:

* approximate k15 threading recovered some otherwise-unthreaded reads but was not
  an exact-threading superset, so its path resolver could lose evidence that was
  already available at k31;
* seed-local rare rescue was directionally useful, but the two-round recruiter
  expanded from only 664 to 685 read pairs, leaving most mate-borne locus flanks
  unused.

Stage17 changes mechanisms, not metric gates:

1. exact-preserving hierarchical threading. Every read first uses exact k31
   threading. Only reads with no exact segment fall back to sparse k19 chaining,
   then k15. The resulting raw-context set is therefore an exact superset. New
   graph paths are add-only on top of Stage10 strict and must be supported by a
   fallback context plus physical graph edges.
2. iterative locus frontier recruitment. Trusted Stage10 rare seeds make the
   initial assignments. Non-backbone k19/k17 kmers from assigned read pairs,
   including their mates, become locus-specific frontier signatures only when
   they do not occur in another seed locus. Several expansion rounds recruit
   new pairs without enabling global singleton rescue.

No reference is used by this assembly stage. Reference-aware breakpoint analysis
is implemented separately for benchmark diagnosis only.
"""
from __future__ import annotations

import argparse
import json
import math
import resource
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import adaptive_k_local_v2 as ak
import graph_path_phaser as gp
import low_abundance_rescue as lr
import repeat_graph_optimizer as rg
import stage14_amplified_methods as s14
import stage15_structural_recovery as s15
import stage16_root_cause as s16
import stage789_optimizer as s78


@dataclass
class HybridPathEvidence:
    path: tuple[str, ...]
    seq: str
    fresh31: int
    fresh_fraction: float
    fallback_support: int
    fallback_span: int
    physical_edges: int
    total_edges: int


def collect_hybrid_contexts(
    graph: gp.Graph,
    exact_index: gp.KmerIndex,
    fallbacks: list[tuple[str, s16.SparseGraphIndex, int, int]],
    read1: Path,
    read2: Path,
    *,
    max_context: int = 10,
) -> tuple[
    Counter[tuple[str, ...]],
    Counter[tuple[str, ...]],
    Counter[tuple[str, ...]],
    dict[str, object],
]:
    """Thread exact first; approximate chaining is fallback-only.

    The returned hybrid counter is literally exact_ctx + fallback_ctx, so exact
    evidence can never be removed by the approximate method.
    """
    exact_ctx: Counter[tuple[str, ...]] = Counter()
    fallback_ctx: Counter[tuple[str, ...]] = Counter()
    stats: dict[str, object] = {
        "reads": 0,
        "exact_threaded_reads": 0,
        "exact_segments": 0,
        "unthreaded_after_all": 0,
    }
    for label, index, _window, _beam in fallbacks:
        stats[f"{label}_indexed_reads"] = 0
        stats[f"{label}_ambiguous_reads"] = 0
        stats[f"{label}_threaded_reads"] = 0
        stats[f"{label}_graph_bridge_reads"] = 0
        stats[f"{label}_anchors_on_chains"] = 0
        stats[f"{label}_ambiguous_anchors_on_chains"] = 0
        stats[f"{label}_dropped_repetitive_keys"] = index.dropped_repetitive_keys

    caches: dict[str, dict[tuple[str, str], list[str] | None]] = {
        label: {} for label, _index, _window, _beam in fallbacks
    }
    for fastq in (read1, read2):
        for _name, seq in gp.read_fastq(fastq):
            stats["reads"] = int(stats["reads"]) + 1
            exact_segments = gp.thread_sequence(seq, graph, exact_index, None)
            if exact_segments:
                stats["exact_threaded_reads"] = int(stats["exact_threaded_reads"]) + 1
                stats["exact_segments"] = int(stats["exact_segments"]) + len(exact_segments)
                for segment in exact_segments:
                    gp.add_context(exact_ctx, segment, max_context)
                continue

            accepted = False
            for label, index, window, beam in fallbacks:
                events = s16.anchor_events(seq, index, window)
                if events:
                    stats[f"{label}_indexed_reads"] = int(stats[f"{label}_indexed_reads"]) + 1
                ambiguous = sum(len(event.candidates) > 1 for event in events)
                if ambiguous:
                    stats[f"{label}_ambiguous_reads"] = int(stats[f"{label}_ambiguous_reads"]) + 1
                chain = s16.chain_anchor_events(
                    events,
                    graph,
                    beam_width=beam,
                    cache=caches[label],
                )
                if chain is None:
                    continue
                stats[f"{label}_threaded_reads"] = int(stats[f"{label}_threaded_reads"]) + 1
                stats[f"{label}_graph_bridge_reads"] = int(stats[f"{label}_graph_bridge_reads"]) + int(chain.graph_bridges > 0)
                stats[f"{label}_anchors_on_chains"] = int(stats[f"{label}_anchors_on_chains"]) + chain.anchors
                stats[f"{label}_ambiguous_anchors_on_chains"] = int(stats[f"{label}_ambiguous_anchors_on_chains"]) + chain.ambiguous_anchors
                gp.add_context(fallback_ctx, list(chain.path), max_context)
                accepted = True
                break
            if not accepted:
                stats["unthreaded_after_all"] = int(stats["unthreaded_after_all"]) + 1

    hybrid = Counter(exact_ctx)
    hybrid.update(fallback_ctx)
    stats["exact_contexts"] = len(exact_ctx)
    stats["fallback_contexts"] = len(fallback_ctx)
    stats["fallback_context_weight"] = sum(fallback_ctx.values())
    stats["hybrid_contexts"] = len(hybrid)
    stats["hybrid_context_weight"] = sum(hybrid.values())
    for label, _index, _window, _beam in fallbacks:
        stats[f"{label}_bridge_cache_entries"] = len(caches[label])
    return exact_ctx, fallback_ctx, hybrid, stats


def path_context_support(
    path: list[str] | tuple[str, ...],
    contexts: Counter[tuple[str, ...]],
    *,
    max_span: int = 6,
) -> tuple[int, int]:
    best_support = 0
    best_span = 0
    nodes = list(path)
    for span in range(2, min(max_span, len(nodes)) + 1):
        for start in range(0, len(nodes) - span + 1):
            support = contexts.get(tuple(nodes[start : start + span]), 0)
            if span > best_span and support > 0:
                best_support, best_span = support, span
            elif span == best_span and support > best_support:
                best_support = support
    return best_support, best_span


def physical_edge_fraction(graph: gp.Graph, path: list[str] | tuple[str, ...]) -> tuple[int, int]:
    supported = 0
    total = 0
    nodes = list(path)
    for left, right in zip(nodes, nodes[1:]):
        total += 1
        ev = graph.edge.get((left, right), gp.EdgeEvidence())
        supported += int((ev.direct + ev.gapped + ev.pairs) > 0)
    return supported, total


def select_hybrid_path_additions(
    paths: list[list[str]],
    graph: gp.Graph,
    strict_baseline: Path,
    fallback_ctx: Counter[tuple[str, ...]],
) -> list[HybridPathEvidence]:
    baseline31, baseline_bases = lr.backbone_kmers(strict_baseline, 31)
    candidates: list[HybridPathEvidence] = []
    for path in paths:
        if len(path) < 2:
            continue
        support, span = path_context_support(path, fallback_ctx)
        if support <= 0:
            continue
        seq = lr.canonical(gp.path_sequence(path, graph))
        if len(seq) < 200:
            continue
        all31 = set(lr.kmers(seq, 31))
        fresh = all31 - baseline31
        if len(fresh) < 8:
            continue
        physical, total = physical_edge_fraction(graph, path)
        if total and physical / total < 0.80:
            continue
        candidates.append(
            HybridPathEvidence(
                path=tuple(path),
                seq=seq,
                fresh31=len(fresh),
                fresh_fraction=len(fresh) / max(1, len(all31)),
                fallback_support=support,
                fallback_span=span,
                physical_edges=physical,
                total_edges=total,
            )
        )

    candidates.sort(
        key=lambda item: (
            -item.fallback_span,
            -item.fallback_support,
            -item.fresh31,
            -item.fresh_fraction,
            -len(item.seq),
            item.seq,
        )
    )
    selected: list[HybridPathEvidence] = []
    selected_fresh: set[str] = set()
    total_bases = 0
    max_bases = max(100_000, int(baseline_bases * 0.08))
    for item in candidates:
        fresh = (set(lr.kmers(item.seq, 31)) - baseline31) - selected_fresh
        if len(fresh) < max(8, math.ceil(0.25 * item.fresh31)):
            continue
        if total_bases + len(item.seq) > max_bases:
            continue
        selected.append(item)
        selected_fresh.update(fresh)
        total_bases += len(item.seq)
    return selected


def write_hybrid_evidence(items: list[HybridPathEvidence], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "name\tlength\tfresh31\tfresh_fraction\tfallback_support\t"
            "fallback_span\tphysical_edges\ttotal_edges\tpath\n"
        )
        for idx, item in enumerate(items, 1):
            handle.write(
                f"hybrid_{idx:06d}\t{len(item.seq)}\t{item.fresh31}\t"
                f"{item.fresh_fraction:.6f}\t{item.fallback_support}\t"
                f"{item.fallback_span}\t{item.physical_edges}\t{item.total_edges}\t"
                f"{','.join(item.path)}\n"
            )


def build_hybrid_thread_candidate(
    scripts: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    read1: Path,
    read2: Path,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    base = pipeline_dir / "current_pipeline"
    graph_opt = pipeline_dir / "graph_optimizer"
    repeat_opt = pipeline_dir / "repeat_optimizer"
    target_gfa = base / "iterative" / "k31_resolve" / "assembly.gfa"
    projection_primary = base / "iterative" / "k21_recall" / "primary_contigs.fasta"
    projection_haplotigs = base / "iterative" / "k21_recall" / "haplotigs.fasta"
    highk_gfa = base / "iterative" / "k55_resolve" / "assembly.gfa"
    base_paths = graph_opt / "stage4_second_pass.paths.tsv"

    graph = gp.Graph.from_gfa(target_gfa)
    membership = gp.preliminary_membership(rg.load_paths(base_paths))
    exact_index = gp.KmerIndex(graph, 31)
    fallback19 = s16.SparseGraphIndex(graph, 19, max_occurrences=5)
    fallback15 = s16.SparseGraphIndex(graph, 15, max_occurrences=4)

    started = time.monotonic()
    exact_ctx, fallback_ctx, hybrid_ctx, hybrid_stats = collect_hybrid_contexts(
        graph,
        exact_index,
        [
            ("k19", fallback19, 6, 8),
            ("k15", fallback15, 8, 6),
        ],
        read1,
        read2,
        max_context=10,
    )
    timings["hybrid_thread_full_library"] = time.monotonic() - started

    proj_ctx, high_ctx, projection_stats = rg.collect_projection_contexts(
        graph,
        exact_index,
        [projection_primary, projection_haplotigs],
        [highk_gfa],
        repeat_opt,
        8,
    )
    second_ctx, second_stats = gp.collect_read_contexts(
        graph, exact_index, read1, read2, membership, 8
    )
    for key in list(second_ctx):
        base_support = exact_ctx.get(key, 0)
        if second_ctx[key] <= base_support:
            del second_ctx[key]
        else:
            second_ctx[key] -= base_support
    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph, exact_index, read1, read2, membership, 8, 8, 420
    )
    repeat_ctx = rg.combined_contexts(second_ctx, pair_ctx)
    all_ctx = rg.combined_contexts(hybrid_ctx, proj_ctx, high_ctx, repeat_ctx)
    simplified, simplify_stats = rg.simplify_graph(graph, all_ctx)
    paths, resolve_stats = s78.resolve_lookahead_seeded_paths(
        simplified,
        hybrid_ctx,
        proj_ctx,
        high_ctx,
        repeat_ctx,
        0.70,
        4,
        200,
        4,
        5,
        0.70,
        0.58,
        1.10,
    )

    outdir = pipeline_dir / "stage17_frontier_recovery" / "hybrid_thread"
    outdir.mkdir(parents=True, exist_ok=True)
    selected = select_hybrid_path_additions(paths, simplified, strict_baseline, fallback_ctx)
    write_hybrid_evidence(selected, outdir / "hybrid_thread_selected.tsv")
    additions = outdir / "hybrid_thread_additions.fasta"
    s14.write_fasta(
        ((f"stage17_hybrid_{idx:06d}", item.seq) for idx, item in enumerate(selected, 1)),
        additions,
    )
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_hybrid_thread",
        timings,
        min_overlap=81,
    )
    return final, additions, {
        "graph_nodes": len(graph.seqs),
        "graph_edges": len(graph.edge),
        "threading": hybrid_stats,
        "projection": projection_stats,
        "second_pass": second_stats,
        "pair_threading": pair_stats,
        "simplification": simplify_stats,
        "path_resolution": resolve_stats,
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "selected_fresh31": sum(item.fresh31 for item in selected),
    }


def load_read_pairs(read1: Path, read2: Path) -> list[tuple[str, str]]:
    return [
        (left[1], right[1])
        for left, right in zip(ak.fastq_records(read1), ak.fastq_records(read2))
    ]


def initial_pair_assignments(
    pairs: list[tuple[str, str]], signatures: dict[str, int]
) -> tuple[dict[int, int], int]:
    assignments: dict[int, int] = {}
    ambiguous = 0
    for idx, (left, right) in enumerate(pairs):
        sid, rejected = s16.unique_assignment(
            s16.score_pair(left, right, signatures, 21, 2),
            2,
            1,
        )
        ambiguous += int(rejected)
        if sid is not None:
            assignments[idx] = sid
    return assignments, ambiguous


def locus_unique_frontier(
    pairs: list[tuple[str, str]],
    assignments: dict[int, int],
    baseline_mers: set[str],
    *,
    k: int,
    stride: int,
    min_count: int,
) -> dict[str, int]:
    counts_by_seed: dict[int, Counter[str]] = defaultdict(Counter)
    for idx, sid in assignments.items():
        left, right = pairs[idx]
        fragment = set(lr.kmers(left, k, stride)) | set(lr.kmers(right, k, stride))
        for mer in fragment:
            if mer not in baseline_mers:
                counts_by_seed[sid][mer] += 1
    membership: dict[str, set[int]] = defaultdict(set)
    for sid, counts in counts_by_seed.items():
        for mer, count in counts.items():
            if count >= min_count:
                membership[mer].add(sid)
    return {
        mer: next(iter(sids))
        for mer, sids in membership.items()
        if len(sids) == 1
    }


def recruit_frontier_round(
    pairs: list[tuple[str, str]],
    assignments: dict[int, int],
    frontier: dict[str, int],
    *,
    k: int,
    stride: int,
    min_hits: int,
    margin: int,
) -> tuple[int, int]:
    new_assignments: dict[int, int] = {}
    ambiguous = 0
    for idx, (left, right) in enumerate(pairs):
        if idx in assignments:
            continue
        sid, rejected = s16.unique_assignment(
            s16.score_pair(left, right, frontier, k, stride),
            min_hits,
            margin,
        )
        ambiguous += int(rejected)
        if sid is not None:
            new_assignments[idx] = sid
    assignments.update(new_assignments)
    return len(new_assignments), ambiguous


def recruit_iterative_frontier_pairs(
    read1: Path,
    read2: Path,
    seeds: list[tuple[str, str]],
    initial_signatures: dict[str, int],
    baseline: Path,
) -> tuple[dict[int, int], dict[str, object]]:
    pairs = load_read_pairs(read1, read2)
    assignments, initial_ambiguous = initial_pair_assignments(pairs, initial_signatures)
    initial_pairs = len(assignments)
    round_specs = [
        (19, 3, 1, 2, 1),
        (19, 3, 1, 2, 1),
        (17, 3, 2, 3, 1),
        (17, 3, 2, 3, 1),
    ]
    baseline_cache: dict[int, set[str]] = {}
    round_stats: list[dict[str, int]] = []
    for round_id, (k, stride, min_count, min_hits, margin) in enumerate(round_specs, 1):
        if k not in baseline_cache:
            baseline_cache[k], _ = lr.backbone_kmers(baseline, k)
        frontier = locus_unique_frontier(
            pairs,
            assignments,
            baseline_cache[k],
            k=k,
            stride=stride,
            min_count=min_count,
        )
        added, ambiguous = recruit_frontier_round(
            pairs,
            assignments,
            frontier,
            k=k,
            stride=stride,
            min_hits=min_hits,
            margin=margin,
        )
        round_stats.append(
            {
                "round": round_id,
                "k": k,
                "frontier_signatures": len(frontier),
                "new_pairs": added,
                "ambiguous_rejects": ambiguous,
                "total_pairs": len(assignments),
            }
        )

    per_seed = Counter(assignments.values())
    return assignments, {
        "total_library_pairs": len(pairs),
        "seed_loci": len(seeds),
        "initial_pairs": initial_pairs,
        "initial_ambiguous_pair_rejects": initial_ambiguous,
        "frontier_rounds": round_stats,
        "total_recruited_pairs": len(assignments),
        "seeds_with_recruited_pairs": len(per_seed),
        "max_pairs_per_seed": max(per_seed.values(), default=0),
        "median_pairs_per_seed": (
            int(sorted(per_seed.values())[len(per_seed) // 2]) if per_seed else 0
        ),
    }


def build_frontier_local_candidate(
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
    outdir = pipeline_dir / "stage17_frontier_recovery" / "frontier_local"
    outdir.mkdir(parents=True, exist_ok=True)

    seeds, signatures, _signature_sets, signature_stats = s16.build_all_seed_signatures(
        seed_fasta, backbone, k=21, min_signature_kmers=4
    )
    started = time.monotonic()
    pair_to_seed, recruit_stats = recruit_iterative_frontier_pairs(
        read1, read2, seeds, signatures, backbone
    )
    raw1 = outdir / "targeted_raw_R1.fastq.gz"
    raw2 = outdir / "targeted_raw_R2.fastq.gz"
    s16.write_targeted_pool(read1, read2, pair_to_seed, raw1, raw2)
    timings["frontier_recruitment"] = time.monotonic() - started

    pair_counts = Counter(pair_to_seed.values())
    with (outdir / "seed_recruitment.tsv").open("w") as handle:
        handle.write("seed_id\tseed_name\trecruited_pairs\n")
        for sid, (name, _seq) in enumerate(seeds):
            handle.write(f"{sid}\t{name}\t{pair_counts.get(sid, 0)}\n")

    seed_v1 = outdir / "trusted_seed_R1.fastq.gz"
    seed_v2 = outdir / "trusted_seed_R2.fastq.gz"
    timings["frontier_seed_virtualization"] = s15.virtualize(
        scripts,
        [seed_fasta],
        seed_v1,
        seed_v2,
        read_length=91,
        insert_size=190,
        stride=60,
        min_length=190,
    )
    aug1 = outdir / "targeted_aug_R1.fastq.gz"
    aug2 = outdir / "targeted_aug_R2.fastq.gz"
    s15.concat_gzip([raw1, seed_v1], aug1)
    s15.concat_gzip([raw2, seed_v2], aug2)
    raw21_counts = s16.targeted_raw_kmer_support(raw1, raw2, 21)

    assembly_inputs: dict[int, list[Path]] = {}
    profiles: dict[str, object] = {}
    for k in (17, 21, 25, 31):
        asm = outdir / f"k{k}"
        timings[f"frontier_local_k{k}"] = s16.assemble_targeted_k(
            bridgeasm, aug1, aug2, asm, k, threads
        )
        assembly_inputs[k] = [asm / "primary_contigs.fasta", asm / "haplotigs.fasta"]
        profiles[f"k{k}"] = s16.profile(asm / "run_profile.json")

    evidence = s16.local_evidence(
        assembly_inputs, seeds, strict_baseline, raw21_counts
    )
    selected = s16.select_local_evidence(evidence, strict_baseline)
    s16.write_local_evidence(evidence, outdir / "frontier_local_evidence.tsv")
    s16.write_local_evidence(selected, outdir / "frontier_local_selected.tsv")
    additions = outdir / "frontier_local_additions.fasta"
    s14.write_fasta(
        (
            (f"stage17_seed{item.seed_id}_k{item.k}_{idx:06d}", item.seq)
            for idx, item in enumerate(selected, 1)
        ),
        additions,
    )
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        outdir / "candidate_frontier_local",
        timings,
        min_overlap=81,
    )
    selected_seeds = {item.seed_id for item in selected}
    return final, additions, {
        **signature_stats,
        **recruit_stats,
        "recruited_library_fraction": len(pair_to_seed) / max(1, recruit_stats["total_library_pairs"]),
        "raw_supported_k21": len(raw21_counts),
        "evidence_records": len(evidence),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "selected_seed_loci": len(selected_seeds),
        "selected_fresh31": sum(item.fresh31 for item in selected),
        "assembly_profiles": profiles,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    started = time.monotonic()
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
        pipeline / "current_pipeline" / "iterative" / "k31_resolve" / "assembly.gfa",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage17 inputs: " + ", ".join(missing))

    root = pipeline / "stage17_frontier_recovery"
    root.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    hybrid_final, hybrid_additions, hybrid_stats = build_hybrid_thread_candidate(
        scripts,
        pipeline,
        strict_baseline,
        args.read1,
        args.read2,
        timings,
    )
    local_final, local_additions, local_stats = build_frontier_local_candidate(
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
    combined = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [hybrid_additions, local_additions],
        root / "combined" / "candidate_combined",
        timings,
        min_overlap=81,
    )

    stats = {
        "pipeline": "bridge-stage17-frontier-recovery-v1",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "metric_targets": False,
            "hybrid_thread": "exact k31 first, fallback-only k19 then k15; add-only graph-path recovery",
            "frontier_local": "iterative locus-unique mate frontier recruitment; singleton rescue only in targeted pool",
            "sequence_join": "exact overlap >=81 bp",
        },
        "methods": {
            "hybrid_thread": hybrid_stats,
            "frontier_local": local_stats,
        },
        "outputs": {
            "hybrid_thread": str(hybrid_final),
            "frontier_local": str(local_final),
            "combined": str(combined),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (root / "stage17_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
