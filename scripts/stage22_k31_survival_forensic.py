#!/usr/bin/env python3
"""Benchmark-only forensic of where reference-supported k31 paths disappear.

For Stage10 gaps already classified as k21->k31 propagation loss, trace unique
reference k31 nodes and adjacent k31 transitions through:

  raw physical fragments -> BridgeAsm k31 unitig graph -> Stage10 contigs

The assembler remains reference-free. This script is diagnostic only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator

BASE = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name = None
    chunks = []
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


def canonical_key(seq: str) -> int | None:
    fwd = 0
    rev = 0
    for index, ch in enumerate(seq):
        value = BASE.get(ch)
        if value is None:
            return None
        fwd = (fwd << 2) | value
        rev |= (3 - value) << (2 * index)
    return min(fwd, rev)


def rolling_keys(seq: str, k: int) -> Iterator[tuple[int, int]]:
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


def sampled_intervals(start: int, end: int, window: int = 500) -> list[tuple[int, int]]:
    length = end - start
    if length <= 3 * window:
        return [(start, end)]
    mid = start + length // 2
    return [
        (start, start + window),
        (max(start, mid - window // 2), min(end, mid + window // 2)),
        (end - window, end),
    ]


def load_reference_dir(path: Path) -> dict[str, str]:
    refs = {}
    for fasta in sorted(path.iterdir()):
        if fasta.suffix.lower() not in {".fa", ".fna", ".fasta"}:
            continue
        for name, seq in fasta_records(fasta):
            refs[name] = seq
    return refs


def load_target_paths(gaps_tsv: Path, refs: dict[str, str], species: str) -> list[dict]:
    paths = []
    with gaps_tsv.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["species"] != species or row["category"] != "k21_to_k31_propagation_loss":
                continue
            target = row["target"]
            start = int(row["start"])
            end = int(row["end"])
            seq = refs[target]
            for piece_index, (piece_start, piece_end) in enumerate(sampled_intervals(start, end)):
                piece = seq[piece_start:piece_end]
                nodes = [key for _position, key in rolling_keys(piece, 31)]
                edges = [key for _position, key in rolling_keys(piece, 32)]
                if nodes:
                    paths.append(
                        {
                            "target": target,
                            "gap_start": start,
                            "gap_end": end,
                            "piece": piece_index,
                            "piece_start": piece_start,
                            "piece_end": piece_end,
                            "nodes": nodes,
                            "edges": edges,
                        }
                    )
    return paths


def reference_counts(refs: dict[str, str], wanted: set[int], k: int) -> Counter[int]:
    counts: Counter[int] = Counter()
    for seq in refs.values():
        local = {key for _position, key in rolling_keys(seq, k) if key in wanted}
        counts.update(local)
    return counts


def scan_raw(read1: Path, read2: Path, node_targets: set[int], edge_targets: set[int]):
    observations: Counter[int] = Counter()
    fragments: Counter[int] = Counter()
    quality_sum: Counter[int] = Counter()
    edge_observations: Counter[int] = Counter()
    edge_fragments: Counter[int] = Counter()
    pairs = 0
    for left, right in zip(fastq_records(read1), fastq_records(read2)):
        pairs += 1
        fragment_nodes: set[int] = set()
        fragment_edges: set[int] = set()
        for seq, qual in (left, right):
            prefix = [0]
            for ch in qual:
                prefix.append(prefix[-1] + max(0, ord(ch) - 33))
            for position, key in rolling_keys(seq, 31):
                if key in node_targets:
                    observations[key] += 1
                    quality_sum[key] += prefix[position + 31] - prefix[position]
                    fragment_nodes.add(key)
            for _position, key in rolling_keys(seq, 32):
                if key in edge_targets:
                    edge_observations[key] += 1
                    fragment_edges.add(key)
        fragments.update(fragment_nodes)
        edge_fragments.update(fragment_edges)
    return pairs, observations, fragments, quality_sum, edge_observations, edge_fragments


def sequence_seen(paths: list[Path], node_targets: set[int], edge_targets: set[int]):
    nodes: set[int] = set()
    edges: set[int] = set()
    for path in paths:
        for _name, seq in fasta_records(path):
            for _position, key in rolling_keys(seq, 31):
                if key in node_targets:
                    nodes.add(key)
            for _position, key in rolling_keys(seq, 32):
                if key in edge_targets:
                    edges.add(key)
    return nodes, edges


def gfa_seen(path: Path, node_targets: set[int], edge_targets: set[int]):
    nodes: set[int] = set()
    edges: set[int] = set()
    with path.open() as handle:
        for raw in handle:
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] == "*":
                continue
            seq = fields[2]
            for _position, key in rolling_keys(seq, 31):
                if key in node_targets:
                    nodes.add(key)
            for _position, key in rolling_keys(seq, 32):
                if key in edge_targets:
                    edges.add(key)
    return nodes, edges


def percentile(values: list[int], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lo = math.floor(position)
    hi = math.ceil(position)
    if lo == hi:
        return float(ordered[lo])
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--gaps", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--k31-gfa", type=Path, required=True)
    parser.add_argument("--stage10", type=Path, required=True)
    parser.add_argument("--species", default="Pseudomonas_aeruginosa")
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    refs = load_reference_dir(args.reference_dir)
    paths = load_target_paths(args.gaps, refs, args.species)
    all_nodes = {node for path in paths for node in path["nodes"]}
    ref31 = reference_counts(refs, all_nodes, 31)
    unique_nodes = {key for key, count in ref31.items() if count == 1}

    selected_paths = []
    target_nodes: set[int] = set()
    target_edges: set[int] = set()
    for path in paths:
        nodes = path["nodes"]
        edges = path["edges"]
        keep_edges = []
        for index, edge in enumerate(edges):
            keep = index + 1 < len(nodes) and nodes[index] in unique_nodes and nodes[index + 1] in unique_nodes
            keep_edges.append(keep)
            if keep:
                target_nodes.add(nodes[index])
                target_nodes.add(nodes[index + 1])
                target_edges.add(edge)
        selected = dict(path)
        selected["keep_edges"] = keep_edges
        selected_paths.append(selected)

    pairs, observations, fragments, quality_sum, edge_observations, edge_fragments = scan_raw(
        args.read1, args.read2, target_nodes, target_edges
    )
    graph_nodes, graph_edges = gfa_seen(args.k31_gfa, target_nodes, target_edges)
    emitted_nodes, emitted_edges = sequence_seen([args.stage10], target_nodes, target_edges)

    node_buckets = defaultdict(lambda: Counter(total=0, graph=0, emitted=0))
    for key in target_nodes:
        count = observations[key]
        fragment_count = fragments[key]
        mean_quality = quality_sum[key] / (count * 31) if count else 0.0
        if count == 0:
            bucket = "raw0"
        elif fragment_count == 1 and mean_quality >= 30:
            bucket = "fragment1_q30plus"
        elif fragment_count == 1 and mean_quality >= 20:
            bucket = "fragment1_q20_30"
        elif fragment_count == 1:
            bucket = "fragment1_qbelow20"
        elif fragment_count >= 2 and mean_quality >= 20:
            bucket = "fragment2plus_q20plus"
        else:
            bucket = "fragment2plus_qbelow20"
        values = node_buckets[bucket]
        values["total"] += 1
        values["graph"] += int(key in graph_nodes)
        values["emitted"] += int(key in emitted_nodes)

    failures: Counter[str] = Counter()
    edge_support = defaultdict(lambda: Counter(total=0, graph=0, emitted=0))
    per_gap = defaultdict(Counter)
    missing_runs: list[int] = []
    singleton_missing_runs: list[int] = []

    for path in selected_paths:
        nodes = path["nodes"]
        edges = path["edges"]
        keep_edges = path["keep_edges"]
        gap_key = (path["target"], path["gap_start"], path["gap_end"])

        absent_run = 0
        singleton_run = 0
        for index, key in enumerate(nodes):
            relevant = (index > 0 and keep_edges[index - 1]) or (index < len(keep_edges) and keep_edges[index])
            if not relevant:
                continue
            missing = key not in graph_nodes and fragments[key] >= 1
            mean_quality = quality_sum[key] / max(1, observations[key] * 31)
            singleton = missing and fragments[key] == 1 and mean_quality >= 20
            if missing:
                absent_run += 1
            elif absent_run:
                missing_runs.append(absent_run)
                absent_run = 0
            if singleton:
                singleton_run += 1
            elif singleton_run:
                singleton_missing_runs.append(singleton_run)
                singleton_run = 0
        if absent_run:
            missing_runs.append(absent_run)
        if singleton_run:
            singleton_missing_runs.append(singleton_run)

        for index, edge in enumerate(edges):
            if not keep_edges[index]:
                continue
            left = nodes[index]
            right = nodes[index + 1]
            raw_edge_observations = edge_observations[edge]
            raw_edge_fragments = edge_fragments[edge]
            if raw_edge_observations == 0:
                category = "no_raw_edge_observation"
            elif left not in graph_nodes or right not in graph_nodes:
                missing_values = []
                for node in (left, right):
                    if node in graph_nodes:
                        continue
                    count = observations[node]
                    fragment_count = fragments[node]
                    mean_quality = quality_sum[node] / (count * 31) if count else 0.0
                    missing_values.append((fragment_count, mean_quality))
                if any(fragment_count == 1 and mean_quality >= 20 for fragment_count, mean_quality in missing_values):
                    category = "node_filter_singleton"
                elif any(mean_quality < 20 for _fragment_count, mean_quality in missing_values):
                    category = "node_filter_quality"
                else:
                    category = "node_filter_other"
            elif edge not in graph_edges:
                category = "graph_edge_missing"
            elif edge not in emitted_edges:
                category = "post_graph_path_loss"
            else:
                category = "survived_stage10"
            failures[category] += 1
            per_gap[gap_key][category] += 1
            per_gap[gap_key]["total"] += 1

            support_bucket = (
                "raw_edge_1frag"
                if raw_edge_fragments == 1
                else "raw_edge_2plus_frag"
                if raw_edge_fragments >= 2
                else "raw_edge_0frag"
            )
            values = edge_support[support_bucket]
            values["total"] += 1
            values["graph"] += int(edge in graph_edges)
            values["emitted"] += int(edge in emitted_edges)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = prefix.with_name(prefix.name + "_summary.tsv")
    with summary_path.open("w") as handle:
        handle.write("section\tcategory\ttotal\tgraph_present\tgraph_fraction\temitted\temitted_fraction\n")
        for category, values in sorted(node_buckets.items()):
            total = values["total"]
            handle.write(
                f"node_support\t{category}\t{total}\t{values['graph']}\t{values['graph']/max(1,total):.6f}\t"
                f"{values['emitted']}\t{values['emitted']/max(1,total):.6f}\n"
            )
        for category, values in sorted(edge_support.items()):
            total = values["total"]
            handle.write(
                f"edge_support\t{category}\t{total}\t{values['graph']}\t{values['graph']/max(1,total):.6f}\t"
                f"{values['emitted']}\t{values['emitted']/max(1,total):.6f}\n"
            )
        total_failures = sum(failures.values())
        for category, count in failures.most_common():
            handle.write(
                f"first_failure\t{category}\t{count}\t0\t{count/max(1,total_failures):.6f}\t0\t0.000000\n"
            )
        for name, values in (
            ("raw_observed_graph_missing_run", missing_runs),
            ("highq_singleton_graph_missing_run", singleton_missing_runs),
        ):
            handle.write(f"run_length\t{name}:runs\t{len(values)}\t0\t0\t0\t0\n")
            handle.write(f"run_length\t{name}:median\t{percentile(values,0.5):.3f}\t0\t0\t0\t0\n")
            handle.write(f"run_length\t{name}:p90\t{percentile(values,0.9):.3f}\t0\t0\t0\t0\n")
            handle.write(f"run_length\t{name}:max\t{max(values) if values else 0}\t0\t0\t0\t0\n")

    gap_path = prefix.with_name(prefix.name + "_gaps.tsv")
    categories = [
        "no_raw_edge_observation",
        "node_filter_singleton",
        "node_filter_quality",
        "node_filter_other",
        "graph_edge_missing",
        "post_graph_path_loss",
        "survived_stage10",
    ]
    with gap_path.open("w") as handle:
        handle.write("target\tstart\tend\ttotal\t" + "\t".join(categories) + "\n")
        for (target, start, end), counts in sorted(per_gap.items()):
            handle.write(
                f"{target}\t{start}\t{end}\t{counts['total']}\t"
                + "\t".join(str(counts[category]) for category in categories)
                + "\n"
            )

    stats = {
        "species": args.species,
        "read_pairs": pairs,
        "sampled_pieces": len(selected_paths),
        "unique_target_nodes": len(target_nodes),
        "unique_target_edges": len(target_edges),
        "node_support": {key: dict(value) for key, value in node_buckets.items()},
        "edge_support": {key: dict(value) for key, value in edge_support.items()},
        "first_failure": dict(failures),
        "first_failure_fraction": {
            key: value / max(1, sum(failures.values())) for key, value in failures.items()
        },
        "raw_observed_graph_missing_runs": {
            "count": len(missing_runs),
            "median": percentile(missing_runs, 0.5),
            "p90": percentile(missing_runs, 0.9),
            "max": max(missing_runs) if missing_runs else 0,
        },
        "highq_singleton_graph_missing_runs": {
            "count": len(singleton_missing_runs),
            "median": percentile(singleton_missing_runs, 0.5),
            "p90": percentile(singleton_missing_runs, 0.9),
            "max": max(singleton_missing_runs) if singleton_missing_runs else 0,
        },
    }
    json_path = prefix.with_name(prefix.name + "_stats.json")
    json_path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(summary_path.read_text(), end="")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
