#!/usr/bin/env python3
"""Stage24 benchmark-only forensic for no-k32 Pseudomonas gaps.

Stage22 showed that most target k31 transitions in propagation-loss regions have
no exact raw k32 observation.  Such transitions cannot be repaired by relaxing
k31 node filtering.  This diagnostic asks a narrower question before changing
the assembler:

  When a consecutive run of target k32 transitions is absent from the raw
  library, is the same spelling still present as a continuous path in the
  existing k21 graph?

For each consecutive no-raw-k32 run we reconstruct the directed k21 de Bruijn
adjacency represented by the k21 GFA and classify it using only benchmark truth:

* exact_nonbranching_both_k31_anchors: entire spelling exists at k21, the
  internal k21 path is non-branching, and both boundary k31 nodes survive.
* exact_branching_both_k31_anchors: entire spelling exists, but lower-k graph
  ambiguity must be resolved by read/pair context.
* exact_missing_k31_anchor: lower-k spelling exists but one/both k31 anchors
  are missing, so direct cross-k bridging is not yet possible.
* k21_graph_broken: most k21 nodes exist but one or more required adjacencies
  are absent.
* k21_sequence_missing: substantial lower-k sequence is absent too.

The reference is used only after assembly to measure recoverability.  No output
from this script is consumed by production assembly decisions.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator

import stage22_k31_survival_forensic as s22

BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "a": 0,
    "c": 1,
    "g": 2,
    "t": 3,
}


def oriented_keys(seq: str, k: int) -> Iterator[int]:
    """Yield forward-orientation 2-bit k-mers, resetting across invalid bases."""
    mask = (1 << (2 * k)) - 1
    key = 0
    valid = 0
    for ch in seq:
        value = BASE.get(ch)
        if value is None:
            key = 0
            valid = 0
            continue
        key = ((key << 2) | value) & mask
        valid += 1
        if valid >= k:
            yield key


def graph_topology(path: Path, k: int = 21):
    """Reconstruct directed k-mer adjacency from compacted GFA S records.

    BridgeAsm unitigs spell every raw DBG transition in an S record; L links
    overlap by k bases and therefore do not introduce an additional (k+1)-mer
    that is absent from the unitig spellings.
    """
    nodes: set[int] = set()
    edges: set[tuple[int, int]] = set()
    out: dict[int, set[int]] = defaultdict(set)
    inc: dict[int, set[int]] = defaultdict(set)
    sequences = 0
    with path.open() as handle:
        for raw in handle:
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] == "*":
                continue
            sequences += 1
            keys = list(oriented_keys(fields[2], k))
            nodes.update(keys)
            for left, right in zip(keys, keys[1:]):
                edges.add((left, right))
                out[left].add(right)
                inc[right].add(left)
    return nodes, edges, out, inc, sequences


def consecutive_true_runs(mask: list[bool]) -> list[tuple[int, int]]:
    """Return inclusive [start,end] runs of True values."""
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, value in enumerate(mask):
        if value and start is None:
            start = index
        elif not value and start is not None:
            runs.append((start, index - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


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


def evaluate_k21_spelling(
    seq: str,
    nodes: set[int],
    edges: set[tuple[int, int]],
    out: dict[int, set[int]],
    inc: dict[int, set[int]],
    *,
    k: int = 21,
) -> dict[str, int | float | bool]:
    keys = list(oriented_keys(seq, k))
    transitions = list(zip(keys, keys[1:]))
    node_present = sum(key in nodes for key in keys)
    edge_present = sum(edge in edges for edge in transitions)
    internal = keys[1:-1]
    branch_nodes = sum(
        len(out.get(key, set())) > 1 or len(inc.get(key, set())) > 1
        for key in internal
    )
    return {
        "k21_nodes": len(keys),
        "k21_nodes_present": node_present,
        "k21_node_fraction": node_present / max(1, len(keys)),
        "k21_transitions": len(transitions),
        "k21_transitions_present": edge_present,
        "k21_transition_fraction": edge_present / max(1, len(transitions)),
        "k21_exact_path": bool(transitions) and edge_present == len(transitions),
        "k21_internal_branch_nodes": branch_nodes,
        "k21_nonbranching": bool(transitions) and edge_present == len(transitions) and branch_nodes == 0,
    }


def classify_run(metrics: dict[str, int | float | bool], both_anchors: bool) -> str:
    edge_fraction = float(metrics["k21_transition_fraction"])
    node_fraction = float(metrics["k21_node_fraction"])
    exact = bool(metrics["k21_exact_path"])
    nonbranching = bool(metrics["k21_nonbranching"])
    if exact and both_anchors and nonbranching:
        return "exact_nonbranching_both_k31_anchors"
    if exact and both_anchors:
        return "exact_branching_both_k31_anchors"
    if exact:
        return "exact_missing_k31_anchor"
    if node_fraction >= 0.80 or edge_fraction >= 0.80:
        return "k21_graph_broken"
    return "k21_sequence_missing"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-dir", type=Path, required=True)
    ap.add_argument("--gaps", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--k21-gfa", type=Path, required=True)
    ap.add_argument("--k31-gfa", type=Path, required=True)
    ap.add_argument("--species", default="Pseudomonas_aeruginosa")
    ap.add_argument("--output-prefix", type=Path, required=True)
    args = ap.parse_args()

    refs = s22.load_reference_dir(args.reference_dir)
    paths = s22.load_target_paths(args.gaps, refs, args.species)
    all_nodes31 = {node for path in paths for node in path["nodes"]}
    ref31 = s22.reference_counts(refs, all_nodes31, 31)
    unique31 = {key for key, count in ref31.items() if count == 1}

    selected_paths: list[dict[str, object]] = []
    target_nodes31: set[int] = set()
    target_edges32: set[int] = set()
    for path in paths:
        nodes31 = path["nodes"]
        edges32 = path["edges"]
        keep_edges: list[bool] = []
        for index, edge in enumerate(edges32):
            keep = (
                index + 1 < len(nodes31)
                and nodes31[index] in unique31
                and nodes31[index + 1] in unique31
            )
            keep_edges.append(keep)
            if keep:
                target_nodes31.add(nodes31[index])
                target_nodes31.add(nodes31[index + 1])
                target_edges32.add(edge)
        item = dict(path)
        item["keep_edges"] = keep_edges
        selected_paths.append(item)

    (
        read_pairs,
        _node_observations,
        _node_fragments,
        _node_quality,
        edge_observations,
        edge_fragments,
    ) = s22.scan_raw(args.read1, args.read2, target_nodes31, target_edges32)
    graph31_nodes, _ = s22.gfa_seen(args.k31_gfa, target_nodes31, set())
    k21_nodes, k21_edges, k21_out, k21_inc, k21_s_records = graph_topology(args.k21_gfa, 21)

    rows: list[dict[str, object]] = []
    categories: Counter[str] = Counter()
    category_runs: Counter[str] = Counter()
    exact_lengths: list[int] = []
    actionable_lengths: list[int] = []
    no_raw_transitions = 0
    one_fragment_edges = 0
    two_plus_fragment_edges = 0

    for path in selected_paths:
        nodes31 = path["nodes"]
        edges32 = path["edges"]
        keep_edges = path["keep_edges"]
        no_raw_mask = [
            bool(keep_edges[index]) and edge_observations[edge] == 0
            for index, edge in enumerate(edges32)
        ]
        for index, edge in enumerate(edges32):
            if not keep_edges[index] or edge_observations[edge] == 0:
                continue
            if edge_fragments[edge] == 1:
                one_fragment_edges += 1
            elif edge_fragments[edge] >= 2:
                two_plus_fragment_edges += 1
        piece_seq = refs[str(path["target"])][int(path["piece_start"]):int(path["piece_end"])]
        for run_index, (start, end) in enumerate(consecutive_true_runs(no_raw_mask)):
            transition_count = end - start + 1
            no_raw_transitions += transition_count
            # A transition at i is the 32-mer piece[i:i+32].  Consecutive
            # transitions i..j therefore spell piece[i:j+32].
            seq = piece_seq[start : end + 32]
            left_anchor = nodes31[start] in graph31_nodes
            right_anchor = nodes31[end + 1] in graph31_nodes
            metrics = evaluate_k21_spelling(seq, k21_nodes, k21_edges, k21_out, k21_inc)
            category = classify_run(metrics, left_anchor and right_anchor)
            categories[category] += transition_count
            category_runs[category] += 1
            if bool(metrics["k21_exact_path"]):
                exact_lengths.append(transition_count)
            if category == "exact_nonbranching_both_k31_anchors":
                actionable_lengths.append(transition_count)
            rows.append(
                {
                    "target": path["target"],
                    "gap_start": path["gap_start"],
                    "gap_end": path["gap_end"],
                    "piece": path["piece"],
                    "piece_start": path["piece_start"],
                    "piece_end": path["piece_end"],
                    "run_index": run_index,
                    "transition_start": start,
                    "transition_end": end,
                    "no_raw_k32_transitions": transition_count,
                    "spelling_bases": len(seq),
                    "left_k31_anchor": int(left_anchor),
                    "right_k31_anchor": int(right_anchor),
                    **metrics,
                    "category": category,
                }
            )

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    runs_path = prefix.with_name(prefix.name + "_runs.tsv")
    summary_path = prefix.with_name(prefix.name + "_summary.tsv")
    stats_path = prefix.with_name(prefix.name + "_stats.json")

    fields = list(rows[0]) if rows else ["category"]
    with runs_path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in ("k21_node_fraction", "k21_transition_fraction"):
                if key in out:
                    out[key] = f"{float(out[key]):.6f}"
            writer.writerow(out)

    with summary_path.open("w") as handle:
        handle.write("category\truns\tno_raw_k32_transitions\ttransition_fraction\n")
        for category, count in categories.most_common():
            handle.write(
                f"{category}\t{category_runs[category]}\t{count}\t{count/max(1,no_raw_transitions):.6f}\n"
            )

    stats = {
        "species": args.species,
        "read_pairs": read_pairs,
        "sampled_paths": len(selected_paths),
        "target_unique_k31_nodes": len(target_nodes31),
        "target_unique_k32_transitions": len(target_edges32),
        "no_raw_k32_transitions": no_raw_transitions,
        "raw_observed_k32_one_fragment": one_fragment_edges,
        "raw_observed_k32_two_plus_fragments": two_plus_fragment_edges,
        "no_raw_runs": len(rows),
        "k21_graph_s_records": k21_s_records,
        "k21_graph_oriented_nodes": len(k21_nodes),
        "k21_graph_oriented_edges": len(k21_edges),
        "category_transitions": dict(categories),
        "category_runs": dict(category_runs),
        "exact_path_transition_fraction": (
            sum(categories[name] for name in categories if name.startswith("exact_"))
            / max(1, no_raw_transitions)
        ),
        "both_anchor_exact_transition_fraction": (
            (categories["exact_nonbranching_both_k31_anchors"] + categories["exact_branching_both_k31_anchors"])
            / max(1, no_raw_transitions)
        ),
        "direct_actionable_transition_fraction": (
            categories["exact_nonbranching_both_k31_anchors"] / max(1, no_raw_transitions)
        ),
        "exact_run_lengths_transitions": {
            "median": percentile(exact_lengths, 0.50),
            "p90": percentile(exact_lengths, 0.90),
            "max": max(exact_lengths, default=0),
        },
        "actionable_run_lengths_transitions": {
            "median": percentile(actionable_lengths, 0.50),
            "p90": percentile(actionable_lengths, 0.90),
            "max": max(actionable_lengths, default=0),
        },
        "policy": {
            "benchmark_reference_only": True,
            "production_feedback": False,
            "direct_actionable_definition": "entire no-k32 run is an exact nonbranching k21 path with both boundary k31 nodes present",
        },
    }
    stats_path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(summary_path.read_text())
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
