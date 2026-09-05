#!/usr/bin/env python3
"""Stage23: rescue the k31 singleton-node failure identified by Stage22.

Stage22 showed that, conditional on a reference transition actually occurring in
raw reads, ~80% of the Pseudomonas transitions first disappear because one of
the k31 nodes has only one physical-fragment observation.  Global singleton
rescue is unsafe (Stage15), so this experiment makes the rescue locus-specific
and topology-conservative.

A target-k singleton is eligible only when:
  * it occurs in exactly one real physical read pair and has mean Q >= 30;
  * it lies on an internally non-branching previous-k unitig;
  * the selected previous-k segment is bounded by >=2 target-k solid nodes;
  * every target-(k+1) transition in the selected segment was observed in a
    real physical read pair; and
  * the singleton nodes in the segment are supported by >=2 distinct physical
    read pairs, preventing one erroneous fragment from creating a path.

All selected segments are concatenated with N^k separators into ONE synthetic
physical fragment at Q20.  Therefore any target-k k-mer receives at most +1
fragment support globally.  Because every selected edge was already observed
in raw reads, this prior can repair the node-count gate but cannot invent graph
topology.

The assembler itself remains reference-free.  References are used only by the
benchmark workflow to evaluate this candidate and compare breakpoint retention.
"""
from __future__ import annotations

import argparse
import gzip
import json
import resource
import shutil
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import low_abundance_rescue as lr
import stage14_amplified_methods as s14
import stage18_short_handoff as s18

BASE = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}


@dataclass
class PriorRecord:
    name: str
    seq: str


@dataclass
class GuidedSegment:
    name: str
    seq: str
    singleton_keys: set[int]
    singleton_fragments: set[int]
    solid_nodes: int


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def fastq_records(path: Path) -> Iterator[tuple[str, str]]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            if not qual:
                raise ValueError(f"truncated FASTQ: {path}")
            yield seq, qual


def rolling_keys(seq: str, k: int) -> Iterator[tuple[int, int]]:
    if k <= 0:
        return
    mask = (1 << (2 * k)) - 1
    high = 2 * (k - 1)
    fwd = rev = valid = 0
    for index, ch in enumerate(seq):
        value = BASE.get(ch)
        if value is None:
            fwd = rev = valid = 0
            continue
        fwd = ((fwd << 2) | value) & mask
        rev = (rev >> 2) | ((3 - value) << high)
        valid += 1
        if valid >= k:
            yield index - k + 1, min(fwd, rev)


def ordered_keys(seq: str, k: int) -> list[int]:
    observed = list(rolling_keys(seq, k))
    expected = max(0, len(seq) - k + 1)
    if len(observed) != expected or any(pos != index for index, (pos, _key) in enumerate(observed)):
        return []
    return [key for _position, key in observed]


def canonical_seq(seq: str) -> str:
    seq = seq.upper()
    rc = seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    return min(seq, rc)


def collect_prior_targets(prior: Path, k: int) -> tuple[list[PriorRecord], set[int], set[int]]:
    records: list[PriorRecord] = []
    nodes: set[int] = set()
    edges: set[int] = set()
    seen: set[str] = set()
    for name, seq0 in fasta_records(prior):
        seq = seq0.upper()
        if len(seq) < k + 1 or "N" in seq:
            continue
        canonical = canonical_seq(seq)
        if canonical in seen:
            continue
        seen.add(canonical)
        node_list = ordered_keys(seq, k)
        edge_list = ordered_keys(seq, k + 1)
        if not node_list or len(edge_list) + 1 != len(node_list):
            continue
        records.append(PriorRecord(name, seq))
        nodes.update(node_list)
        edges.update(edge_list)
    return records, nodes, edges


def scan_raw_support(
    read1: Path,
    read2: Path,
    node_targets: set[int],
    edge_targets: set[int],
    k: int,
) -> tuple[Counter[int], Counter[int], Counter[int], dict[int, int], Counter[int], int]:
    observations: Counter[int] = Counter()
    fragment_count: Counter[int] = Counter()
    quality_sum: Counter[int] = Counter()
    singleton_supporter: dict[int, int] = {}
    edge_fragments: Counter[int] = Counter()
    pair_count = 0

    for left, right in zip(fastq_records(read1), fastq_records(read2)):
        pair_count += 1
        fragment_nodes: set[int] = set()
        fragment_edges: set[int] = set()
        for seq, qual in (left, right):
            prefix = [0]
            for ch in qual:
                prefix.append(prefix[-1] + max(0, ord(ch) - 33))
            for pos, key in rolling_keys(seq, k):
                if key not in node_targets:
                    continue
                observations[key] += 1
                quality_sum[key] += prefix[pos + k] - prefix[pos]
                fragment_nodes.add(key)
            for _pos, key in rolling_keys(seq, k + 1):
                if key in edge_targets:
                    fragment_edges.add(key)

        for key in fragment_nodes:
            before = fragment_count[key]
            fragment_count[key] = before + 1
            if before == 0:
                singleton_supporter[key] = pair_count
            elif before == 1:
                singleton_supporter[key] = -1
        edge_fragments.update(fragment_edges)

    return observations, fragment_count, quality_sum, singleton_supporter, edge_fragments, pair_count


def mean_quality(key: int, observations: Counter[int], quality_sum: Counter[int], k: int) -> float:
    count = observations.get(key, 0)
    if count <= 0:
        return 0.0
    return quality_sum.get(key, 0) / (count * k)


def node_state(
    key: int,
    observations: Counter[int],
    fragment_count: Counter[int],
    quality_sum: Counter[int],
    k: int,
    singleton_quality: float,
) -> int:
    fragments = fragment_count.get(key, 0)
    quality = mean_quality(key, observations, quality_sum, k)
    if fragments >= 2 and quality >= 20.0:
        return 2  # production solid
    if fragments == 1 and quality >= singleton_quality:
        return 1  # Stage22 high-quality singleton failure class
    return 0


def select_guided_segments(
    records: list[PriorRecord],
    *,
    k: int,
    observations: Counter[int],
    fragment_count: Counter[int],
    quality_sum: Counter[int],
    singleton_supporter: dict[int, int],
    edge_fragments: Counter[int],
    singleton_quality: float = 30.0,
    min_distinct_singleton_fragments: int = 2,
) -> tuple[list[GuidedSegment], set[int], dict[str, int | float]]:
    selected: list[GuidedSegment] = []
    rescue_keys: set[int] = set()
    seen_sequences: set[str] = set()
    candidate_components = 0
    rejected_no_two_anchors = 0
    rejected_single_fragment = 0
    rejected_no_singleton = 0

    for record in records:
        nodes = ordered_keys(record.seq, k)
        edges = ordered_keys(record.seq, k + 1)
        if not nodes or len(edges) + 1 != len(nodes):
            continue
        states = [
            node_state(
                key,
                observations,
                fragment_count,
                quality_sum,
                k,
                singleton_quality,
            )
            for key in nodes
        ]

        index = 0
        while index < len(nodes):
            if states[index] == 0:
                index += 1
                continue
            start = index
            while (
                index + 1 < len(nodes)
                and states[index + 1] != 0
                and edge_fragments.get(edges[index], 0) >= 1
            ):
                index += 1
            end = index
            candidate_components += 1

            solids = [pos for pos in range(start, end + 1) if states[pos] == 2]
            if len(solids) < 2:
                rejected_no_two_anchors += 1
                index += 1
                continue
            left = solids[0]
            right = solids[-1]
            singleton_positions = [pos for pos in range(left, right + 1) if states[pos] == 1]
            if not singleton_positions:
                rejected_no_singleton += 1
                index += 1
                continue
            supporters = {
                singleton_supporter.get(nodes[pos], -1)
                for pos in singleton_positions
                if singleton_supporter.get(nodes[pos], -1) > 0
            }
            if len(supporters) < min_distinct_singleton_fragments:
                rejected_single_fragment += 1
                index += 1
                continue

            seq = record.seq[left : right + k]
            canonical = canonical_seq(seq)
            if canonical in seen_sequences:
                index += 1
                continue
            seen_sequences.add(canonical)
            singleton_keys = {nodes[pos] for pos in singleton_positions}
            rescue_keys.update(singleton_keys)
            selected.append(
                GuidedSegment(
                    name=f"{record.name}:{left}-{right+k}",
                    seq=seq,
                    singleton_keys=singleton_keys,
                    singleton_fragments=supporters,
                    solid_nodes=len(solids),
                )
            )
            index += 1

    lengths = sorted(len(item.seq) for item in selected)
    stats: dict[str, int | float] = {
        "prior_records": len(records),
        "candidate_components": candidate_components,
        "selected_segments": len(selected),
        "selected_bases": sum(lengths),
        "selected_singleton_keys": len(rescue_keys),
        "selected_singleton_physical_fragments": len(
            {frag for item in selected for frag in item.singleton_fragments}
        ),
        "rejected_no_two_anchors": rejected_no_two_anchors,
        "rejected_single_fragment": rejected_single_fragment,
        "rejected_no_singleton": rejected_no_singleton,
        "median_segment_length": float(lengths[len(lengths) // 2]) if lengths else 0.0,
        "max_segment_length": max(lengths, default=0),
    }
    return selected, rescue_keys, stats


def write_guided_segments(segments: list[GuidedSegment], fasta: Path) -> None:
    fasta.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open("w") as handle:
        for index, item in enumerate(segments, 1):
            handle.write(
                f">guided_{index:07d} source={item.name} singleton_keys={len(item.singleton_keys)} "
                f"singleton_fragments={len(item.singleton_fragments)} solid_nodes={item.solid_nodes}\n"
            )
            for start in range(0, len(item.seq), 80):
                handle.write(item.seq[start : start + 80] + "\n")


def write_single_fragment_prior(
    segments: list[GuidedSegment],
    out1: Path,
    out2: Path,
    *,
    k: int,
    phred: int = 20,
) -> dict[str, int]:
    out1.parent.mkdir(parents=True, exist_ok=True)
    separator = "N" * k
    joined = separator.join(item.seq for item in segments)
    qchar = chr(33 + phred)
    with gzip.open(out1, "wt", compresslevel=3) as left, gzip.open(out2, "wt", compresslevel=3) as right:
        if joined:
            left.write(f"@guided_singleton_prior_k{k}/1\n{joined}\n+\n{qchar * len(joined)}\n")
            right.write(f"@guided_singleton_prior_k{k}/2\nN\n+\n!\n")
    return {
        "synthetic_physical_fragments": int(bool(joined)),
        "synthetic_r1_bases": len(joined),
        "max_synthetic_fragment_support_per_kmer": int(bool(joined)),
        "phred": phred,
    }


def fasta_keyset(path: Path, k: int) -> set[int]:
    result: set[int] = set()
    for _name, seq in fasta_records(path):
        result.update(key for _pos, key in rolling_keys(seq, k))
    return result


def select_guided_additions(
    inputs: list[Path],
    baseline: Path,
    rescue_keys: set[int],
    output: Path,
    *,
    k: int = 31,
) -> dict[str, int]:
    baseline_keys = fasta_keyset(baseline, k)
    seen: set[str] = set()
    selected: list[tuple[str, str, int, int]] = []
    scanned = 0
    for path in inputs:
        for name, seq0 in fasta_records(path):
            scanned += 1
            seq = seq0.upper()
            keys = {key for _pos, key in rolling_keys(seq, k)}
            rescue_hits = len(keys & rescue_keys)
            fresh = len(keys - baseline_keys)
            if rescue_hits == 0 or fresh == 0:
                continue
            canonical = canonical_seq(seq)
            if canonical in seen:
                continue
            seen.add(canonical)
            selected.append((name, seq, rescue_hits, fresh))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for index, (_name, seq, rescue_hits, fresh) in enumerate(
            sorted(selected, key=lambda item: (-len(item[1]), item[1])), 1
        ):
            handle.write(
                f">stage23_guided_{index:07d} len={len(seq)} rescue31={rescue_hits} fresh31={fresh}\n"
            )
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")
    return {
        "scanned_output_records": scanned,
        "selected_output_records": len(selected),
        "selected_output_bases": sum(len(item[1]) for item in selected),
        "selected_rescue31_hits": sum(item[2] for item in selected),
        "selected_fresh31": sum(item[3] for item in selected),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeasm", type=Path, required=True)
    parser.add_argument("--pipeline-dir", type=Path, required=True)
    parser.add_argument("-1", "--read1", type=Path, required=True)
    parser.add_argument("-2", "--read2", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--singleton-quality", type=float, default=30.0)
    parser.add_argument("--min-distinct-singleton-fragments", type=int, default=2)
    args = parser.parse_args()

    pipeline = args.pipeline_dir
    k21_unitigs = pipeline / "current_pipeline" / "iterative" / "k21_recall" / "unitigs.fasta"
    stage10 = pipeline / "stage10_multik_rescue"
    strict_baseline = stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    required = [args.bridgeasm, args.read1, args.read2, k21_unitigs, strict_baseline]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage23 inputs: " + ", ".join(missing))

    root = pipeline / "stage23_guided_singleton"
    root.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    timings: dict[str, float] = {}

    t0 = time.monotonic()
    records, node_targets, edge_targets = collect_prior_targets(k21_unitigs, 31)
    (
        observations,
        fragment_count,
        quality_sum,
        singleton_supporter,
        edge_fragments,
        read_pairs,
    ) = scan_raw_support(args.read1, args.read2, node_targets, edge_targets, 31)
    segments, rescue_keys, selection_stats = select_guided_segments(
        records,
        k=31,
        observations=observations,
        fragment_count=fragment_count,
        quality_sum=quality_sum,
        singleton_supporter=singleton_supporter,
        edge_fragments=edge_fragments,
        singleton_quality=args.singleton_quality,
        min_distinct_singleton_fragments=args.min_distinct_singleton_fragments,
    )
    timings["build_guided_prior"] = time.monotonic() - t0

    segments_fasta = root / "guided_segments.fasta"
    write_guided_segments(segments, segments_fasta)
    prior1 = root / "guided_prior_R1.fastq.gz"
    prior2 = root / "guided_prior_R2.fastq.gz"
    prior_stats = write_single_fragment_prior(segments, prior1, prior2, k=31, phred=20)

    aug1 = root / "aug_R1.fastq.gz"
    aug2 = root / "aug_R2.fastq.gz"
    s18.concat_gzip([args.read1, prior1], aug1)
    s18.concat_gzip([args.read2, prior2], aug2)
    guided_asm = root / "guided_k31"
    timings["guided_k31_assembly"] = s18.assemble_target_k(
        args.bridgeasm, aug1, aug2, guided_asm, 31, args.threads
    )

    additions = root / "guided_additions.fasta"
    addition_stats = select_guided_additions(
        [guided_asm / "primary_contigs.fasta", guided_asm / "haplotigs.fasta"],
        strict_baseline,
        rescue_keys,
        additions,
    )
    scripts = Path(__file__).resolve().parent
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        root / "candidate_guided_singleton",
        timings,
        min_overlap=81,
    )
    shutil.copy2(guided_asm / "primary_contigs.fasta", root / "guided_k31_primary.fasta")

    raw_node_categories = Counter()
    for key in node_targets:
        fragments = fragment_count.get(key, 0)
        quality = mean_quality(key, observations, quality_sum, 31)
        if fragments == 0:
            raw_node_categories["raw0"] += 1
        elif fragments == 1 and quality >= args.singleton_quality:
            raw_node_categories["singleton_q30plus"] += 1
        elif fragments == 1:
            raw_node_categories["singleton_other"] += 1
        elif quality >= 20.0:
            raw_node_categories["solid"] += 1
        else:
            raw_node_categories["multi_low_quality"] += 1

    stats = {
        "pipeline": "bridge-stage23-guided-singleton-v1",
        "policy": {
            "reference_free": True,
            "root_cause": "Stage22 raw-observed transition first failure at singleton node filter",
            "prior": "k21 nonbranching unitigs",
            "singleton_quality": args.singleton_quality,
            "solid_definition": "physical_fragments>=2 and meanQ>=20",
            "edge_requirement": "every selected target-k+1 transition observed in >=1 real physical fragment",
            "anchors": "at least two target-k solid nodes",
            "distinct_singleton_fragments": args.min_distinct_singleton_fragments,
            "synthetic_support": "one global Q20 physical fragment; max +1 fragment per target-k kmer",
            "zero_raw_target_nodes_allowed": False,
            "zero_raw_edges_allowed": False,
        },
        "read_pairs": read_pairs,
        "prior_target_nodes": len(node_targets),
        "prior_target_edges": len(edge_targets),
        "raw_node_categories": dict(raw_node_categories),
        "selection": selection_stats,
        "prior_fragment": prior_stats,
        "additions": addition_stats,
        "guided_k31_profile": s18.profile(guided_asm / "run_profile.json"),
        "outputs": {
            "segments": str(segments_fasta),
            "guided_gfa": str(guided_asm / "assembly.gfa"),
            "guided_primary": str(root / "guided_k31_primary.fasta"),
            "additions": str(additions),
            "final": str(final),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (root / "stage23_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
