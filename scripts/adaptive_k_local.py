#!/usr/bin/env python3
"""Reference-free local adaptive-k experiment for BridgeAsm.

The experiment starts from a low-k BridgeAsm graph, identifies branch-heavy
unitig neighborhoods, routes read pairs back to those neighborhoods with
canonical seed k-mers, and reassembles only the unresolved neighborhoods at a
small set of candidate k values. A higher-k result is promoted only when it
retains most local assembled bases and improves graph simplicity or N50.

This is intentionally an experimental orchestration layer. It does not change
the production assembler and it never uses a reference to choose k or contigs.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass
class Segment:
    name: str
    sequence: str
    coverage: float


@dataclass
class Candidate:
    k: int
    output_dir: Path
    bases: int
    n50: int
    contigs: int
    branches: int


@dataclass
class Neighborhood:
    identifier: int
    nodes: set[str]
    branch_nodes: int
    mean_coverage: float
    pair_count: int = 0


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode)


def fastq_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with open_text(path, "r") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ: {path}")
            yield (
                header.rstrip("\n"),
                sequence.rstrip("\n").upper(),
                plus.rstrip("\n"),
                quality.rstrip("\n"),
            )


def write_record(handle: TextIO, record: tuple[str, str, str, str]) -> None:
    handle.write("\n".join(record) + "\n")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def canonical_kmers(sequence: str, k: int, stride: int = 1) -> Iterator[str]:
    if k <= 0 or len(sequence) < k:
        return
    for start in range(0, len(sequence) - k + 1, stride):
        kmer = sequence[start : start + k]
        if set(kmer) <= {"A", "C", "G", "T"}:
            reverse = reverse_complement(kmer)
            yield min(kmer, reverse)


def parse_gfa(path: Path) -> tuple[dict[str, Segment], dict[str, set[str]], Counter[str], Counter[str]]:
    segments: dict[str, Segment] = {}
    undirected: dict[str, set[str]] = defaultdict(set)
    indegree: Counter[str] = Counter()
    outdegree: Counter[str] = Counter()
    with path.open() as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S" and len(fields) >= 3:
                coverage = 0.0
                for tag in fields[3:]:
                    if tag.startswith("KC:f:"):
                        coverage = float(tag[5:])
                segments[fields[1]] = Segment(fields[1], fields[2].upper(), coverage)
            elif fields[0] == "L" and len(fields) >= 5:
                source = fields[1]
                target = fields[3]
                undirected[source].add(target)
                undirected[target].add(source)
                outdegree[source] += 1
                indegree[target] += 1
    for name in segments:
        undirected.setdefault(name, set())
        indegree.setdefault(name, 0)
        outdegree.setdefault(name, 0)
    return segments, undirected, indegree, outdegree


def expand_neighborhood(seed: str, graph: dict[str, set[str]], radius: int, max_nodes: int) -> set[str]:
    seen = {seed}
    queue: deque[tuple[str, int]] = deque([(seed, 0)])
    while queue and len(seen) < max_nodes:
        node, depth = queue.popleft()
        if depth >= radius:
            continue
        for neighbor in sorted(graph[node]):
            if neighbor in seen:
                continue
            seen.add(neighbor)
            queue.append((neighbor, depth + 1))
            if len(seen) >= max_nodes:
                break
    return seen


def build_neighborhoods(
    segments: dict[str, Segment],
    graph: dict[str, set[str]],
    indegree: Counter[str],
    outdegree: Counter[str],
    radius: int,
    max_nodes: int,
    max_neighborhoods: int,
) -> list[Neighborhood]:
    branch_seeds = [
        name
        for name in segments
        if indegree[name] > 1 or outdegree[name] > 1
    ]
    branch_seeds.sort(
        key=lambda name: (
            -(max(0, indegree[name] - 1) + max(0, outdegree[name] - 1)),
            -len(graph[name]),
            len(segments[name].sequence),
            name,
        )
    )

    groups: list[set[str]] = []
    for seed in branch_seeds:
        expanded = expand_neighborhood(seed, graph, radius, max_nodes)
        overlapping = [index for index, group in enumerate(groups) if group & expanded]
        if overlapping:
            merged = set(expanded)
            for index in reversed(overlapping):
                merged.update(groups.pop(index))
            groups.append(merged)
        else:
            groups.append(expanded)
        if len(groups) >= max_neighborhoods * 2:
            break

    groups.sort(
        key=lambda group: (
            -sum(1 for node in group if indegree[node] > 1 or outdegree[node] > 1),
            -sum(len(segments[node].sequence) for node in group),
        )
    )
    groups = groups[:max_neighborhoods]
    output: list[Neighborhood] = []
    for index, group in enumerate(groups, start=1):
        coverages = [segments[node].coverage for node in group if segments[node].coverage > 0]
        output.append(
            Neighborhood(
                identifier=index,
                nodes=group,
                branch_nodes=sum(1 for node in group if indegree[node] > 1 or outdegree[node] > 1),
                mean_coverage=sum(coverages) / len(coverages) if coverages else 0.0,
            )
        )
    return output


def build_signature_index(
    neighborhoods: list[Neighborhood],
    segments: dict[str, Segment],
    seed_k: int,
    stride: int,
    max_memberships: int,
) -> dict[str, tuple[int, ...]]:
    membership: dict[str, set[int]] = defaultdict(set)
    for neighborhood in neighborhoods:
        local: set[str] = set()
        for node in neighborhood.nodes:
            local.update(canonical_kmers(segments[node].sequence, seed_k, stride))
        for kmer in local:
            membership[kmer].add(neighborhood.identifier)
    return {
        kmer: tuple(sorted(ids))
        for kmer, ids in membership.items()
        if 0 < len(ids) <= max_memberships
    }


def route_pairs(
    read1: Path,
    read2: Path,
    output_dir: Path,
    neighborhoods: list[Neighborhood],
    signature_index: dict[str, tuple[int, ...]],
    seed_k: int,
    stride: int,
    min_hits: int,
    hit_margin: int,
) -> None:
    handles: dict[int, tuple[TextIO, TextIO]] = {}
    by_id = {item.identifier: item for item in neighborhoods}
    for item in neighborhoods:
        directory = output_dir / f"neighborhood_{item.identifier:03d}" / "reads"
        directory.mkdir(parents=True, exist_ok=True)
        handles[item.identifier] = (
            gzip.open(directory / "R1.fastq.gz", "wt"),
            gzip.open(directory / "R2.fastq.gz", "wt"),
        )

    left_iter = fastq_records(read1)
    right_iter = fastq_records(read2)
    pair_number = 0
    while True:
        try:
            left = next(left_iter)
        except StopIteration:
            left = None
        try:
            right = next(right_iter)
        except StopIteration:
            right = None
        if left is None and right is None:
            break
        if left is None or right is None:
            raise ValueError("paired FASTQ files contain different record counts")
        pair_number += 1
        scores: Counter[int] = Counter()
        for sequence in (left[1], right[1]):
            for kmer in canonical_kmers(sequence, seed_k, stride):
                for identifier in signature_index.get(kmer, ()):
                    scores[identifier] += 1
        ranked = scores.most_common(2)
        if not ranked or ranked[0][1] < min_hits:
            continue
        if len(ranked) > 1 and ranked[0][1] - ranked[1][1] < hit_margin:
            continue
        identifier = ranked[0][0]
        out1, out2 = handles[identifier]
        write_record(out1, left)
        write_record(out2, right)
        by_id[identifier].pair_count += 1

    for out1, out2 in handles.values():
        out1.close()
        out2.close()


def gfa_branch_count(path: Path) -> int:
    if not path.exists():
        return 0
    _segments, _graph, indegree, outdegree = parse_gfa(path)
    return sum(1 for node in set(indegree) | set(outdegree) if indegree[node] > 1 or outdegree[node] > 1)


def run_bridge(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    output_dir: Path,
    k: int,
    threads: int,
    min_contig_length: int,
) -> Candidate:
    output_dir.mkdir(parents=True, exist_ok=True)
    mercy_span = 24 if k <= 21 else 16
    command = [
        str(bridgeasm),
        "assemble",
        "-1",
        str(read1),
        "-2",
        str(read2),
        "-o",
        str(output_dir),
        "-k",
        str(k),
        "--min-count",
        "2",
        "--mercy-max-kmers",
        str(mercy_span),
        "--mercy-min-support",
        "1",
        "--mercy-min-quality",
        "25",
        "--min-read-support",
        "2",
        "--min-pair-support",
        "2",
        "--min-primary-support",
        "5",
        "--primary-dominance",
        "0.75",
        "--min-contig-length",
        str(min_contig_length),
        "--threads",
        str(threads),
    ]
    subprocess.run(command, check=True)
    profile = json.loads((output_dir / "run_profile.json").read_text())
    return Candidate(
        k=k,
        output_dir=output_dir,
        bases=int(profile.get("primary_bases", 0)),
        n50=int(profile.get("primary_n50", 0)),
        contigs=int(profile.get("primary_contigs", 0)),
        branches=gfa_branch_count(output_dir / "assembly.gfa"),
    )


def promote_candidate(candidates: list[Candidate], base_k: int, min_base_fraction: float) -> tuple[Candidate, bool]:
    if not candidates:
        raise ValueError("no adaptive-k candidates")
    max_bases = max(item.bases for item in candidates)
    floor = max_bases * min_base_fraction
    eligible = [item for item in candidates if item.bases >= floor]
    best = max(eligible, key=lambda item: (-item.branches, item.n50, item.bases, -item.k))
    base = next((item for item in candidates if item.k == base_k), None)
    if base is None or best.k == base_k:
        return best, False
    improved = best.branches < base.branches or best.n50 >= math.ceil(base.n50 * 1.10)
    return best, improved


def run_helpers(base_fasta: Path, selected: list[Path], output_dir: Path, min_contig_length: int) -> Path:
    scripts = Path(__file__).resolve().parent
    raw = output_dir / "adaptive_exact_union.fasta"
    final = output_dir / "adaptive_contigs.fasta"
    subprocess.run(
        [
            sys.executable,
            str(scripts / "merge_fasta_unique.py"),
            str(raw),
            str(base_fasta),
            *map(str, selected),
            "--min-length",
            str(min_contig_length),
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(scripts / "filter_contained_fasta.py"),
            str(raw),
            str(final),
            "--min-length",
            str(min_contig_length),
            "--stats-json",
            str(output_dir / "containment_stats.json"),
            "--removed-tsv",
            str(output_dir / "contained_removed.tsv"),
        ],
        check=True,
    )
    return final


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeasm", required=True, type=Path)
    parser.add_argument("--read1", required=True, type=Path)
    parser.add_argument("--read2", required=True, type=Path)
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--base-k", type=int, default=21)
    parser.add_argument("--candidate-k", default="25,31,41")
    parser.add_argument("--seed-k", type=int, default=15)
    parser.add_argument("--routing-stride", type=int, default=3)
    parser.add_argument("--signature-stride", type=int, default=3)
    parser.add_argument("--max-signature-memberships", type=int, default=2)
    parser.add_argument("--min-seed-hits", type=int, default=3)
    parser.add_argument("--seed-hit-margin", type=int, default=1)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--max-nodes", type=int, default=80)
    parser.add_argument("--max-neighborhoods", type=int, default=6)
    parser.add_argument("--min-pairs", type=int, default=80)
    parser.add_argument("--min-candidate-base-fraction", type=float, default=0.85)
    parser.add_argument("--min-contig-length", type=int, default=200)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    if not args.bridgeasm.exists():
        raise FileNotFoundError(args.bridgeasm)
    base_gfa = args.base_dir / "assembly.gfa"
    base_fasta = args.base_dir / "primary_contigs.fasta"
    if not base_gfa.exists() or not base_fasta.exists():
        raise FileNotFoundError("base-dir must contain assembly.gfa and primary_contigs.fasta")
    if not 0.0 < args.min_candidate_base_fraction <= 1.0:
        raise ValueError("min-candidate-base-fraction must be in (0, 1]")

    candidate_ks = {args.base_k}
    for value in args.candidate_k.split(","):
        if value.strip():
            candidate_ks.add(int(value))
    candidate_ks = {value for value in candidate_ks if value > 0 and value % 2 == 1}
    if args.base_k not in candidate_ks:
        candidate_ks.add(args.base_k)

    args.output.mkdir(parents=True, exist_ok=True)
    segments, graph, indegree, outdegree = parse_gfa(base_gfa)
    neighborhoods = build_neighborhoods(
        segments,
        graph,
        indegree,
        outdegree,
        args.radius,
        args.max_nodes,
        args.max_neighborhoods,
    )

    if not neighborhoods:
        final = run_helpers(base_fasta, [], args.output, args.min_contig_length)
        (args.output / "summary.json").write_text(
            json.dumps({"neighborhoods": 0, "promoted": 0, "final_fasta": str(final)}, indent=2) + "\n"
        )
        print("no unresolved branch neighborhoods; kept baseline")
        return

    signature_index = build_signature_index(
        neighborhoods,
        segments,
        args.seed_k,
        args.signature_stride,
        args.max_signature_memberships,
    )
    route_pairs(
        args.read1,
        args.read2,
        args.output,
        neighborhoods,
        signature_index,
        args.seed_k,
        args.routing_stride,
        args.min_seed_hits,
        args.seed_hit_margin,
    )

    selected_fastas: list[Path] = []
    choices: list[tuple[int, int, int, int, int, int, bool]] = []
    for item in neighborhoods:
        directory = args.output / f"neighborhood_{item.identifier:03d}"
        if item.pair_count < args.min_pairs:
            choices.append((item.identifier, item.pair_count, args.base_k, 0, 0, 0, False))
            continue
        read1 = directory / "reads" / "R1.fastq.gz"
        read2 = directory / "reads" / "R2.fastq.gz"
        candidates: list[Candidate] = []
        for k in sorted(candidate_ks):
            candidates.append(
                run_bridge(
                    args.bridgeasm,
                    read1,
                    read2,
                    directory / f"k{k}",
                    k,
                    args.threads,
                    args.min_contig_length,
                )
            )
        best, promoted = promote_candidate(candidates, args.base_k, args.min_candidate_base_fraction)
        if promoted:
            selected_fastas.append(best.output_dir / "primary_contigs.fasta")
        choices.append(
            (
                item.identifier,
                item.pair_count,
                best.k,
                best.bases,
                best.n50,
                best.branches,
                promoted,
            )
        )
        with (directory / "candidate_metrics.tsv").open("w") as handle:
            handle.write("k\tprimary_bases\tprimary_n50\tprimary_contigs\tbranch_nodes\n")
            for candidate in candidates:
                handle.write(
                    f"{candidate.k}\t{candidate.bases}\t{candidate.n50}\t{candidate.contigs}\t{candidate.branches}\n"
                )

    final = run_helpers(base_fasta, selected_fastas, args.output, args.min_contig_length)
    with (args.output / "neighborhoods.tsv").open("w") as handle:
        handle.write("neighborhood\tnodes\tbranch_nodes\tmean_coverage\tpairs\n")
        for item in neighborhoods:
            handle.write(
                f"{item.identifier}\t{len(item.nodes)}\t{item.branch_nodes}\t{item.mean_coverage:.6f}\t{item.pair_count}\n"
            )
    with (args.output / "adaptive_k_choices.tsv").open("w") as handle:
        handle.write("neighborhood\tpairs\tchosen_k\tprimary_bases\tprimary_n50\tbranch_nodes\tpromoted\n")
        for row in choices:
            handle.write("\t".join(map(str, row)) + "\n")

    summary = {
        "base_k": args.base_k,
        "candidate_k": sorted(candidate_ks),
        "neighborhoods": len(neighborhoods),
        "routed_pairs": sum(item.pair_count for item in neighborhoods),
        "promoted": sum(1 for row in choices if row[-1]),
        "signature_kmers": len(signature_index),
        "final_fasta": str(final),
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
