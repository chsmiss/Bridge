#!/usr/bin/env python3
"""Second local adaptive-k prototype for BridgeAsm.

V1 deliberately exposed a failure mode: transitively merging overlapping branch
neighborhoods could turn a "local" experiment into two ~1,000-node graph
regions. V2 keeps neighborhoods bounded and non-overlapping. It ranks branch
seeds, performs a capped radius BFS around each seed, skips seeds already
covered by a previous accepted neighborhood, and routes reads only with k-mers
that are specific to one accepted neighborhood.

No reference is used for neighborhood selection, candidate-k scoring, or output
promotion. The script remains an experiment wrapper rather than production
assembler logic.
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
_VALID = frozenset("ACGT")


@dataclass
class Segment:
    name: str
    sequence: str
    coverage: float


@dataclass
class Neighborhood:
    identifier: int
    seed: str
    nodes: set[str]
    branch_nodes: int
    mean_coverage: float
    pair_count: int = 0


@dataclass
class Candidate:
    k: int
    output_dir: Path
    bases: int
    n50: int
    contigs: int
    branches: int


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
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record in {path}")
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


def canonical_kmers(sequence: str, k: int, stride: int) -> Iterator[str]:
    if k <= 0 or stride <= 0 or len(sequence) < k:
        return
    for start in range(0, len(sequence) - k + 1, stride):
        kmer = sequence[start : start + k]
        if set(kmer) <= _VALID:
            reverse = reverse_complement(kmer)
            yield min(kmer, reverse)


def parse_gfa(
    path: Path,
) -> tuple[dict[str, Segment], dict[str, set[str]], Counter[str], Counter[str]]:
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


def branch_excess(name: str, indegree: Counter[str], outdegree: Counter[str]) -> int:
    return max(0, indegree[name] - 1) + max(0, outdegree[name] - 1)


def expand_neighborhood(
    seed: str,
    graph: dict[str, set[str]],
    radius: int,
    max_nodes: int,
) -> set[str]:
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
    min_branch_nodes: int,
    max_overlap_fraction: float,
) -> list[Neighborhood]:
    seeds = [name for name in segments if branch_excess(name, indegree, outdegree) > 0]
    seeds.sort(
        key=lambda name: (
            -branch_excess(name, indegree, outdegree),
            -len(graph[name]),
            -segments[name].coverage,
            len(segments[name].sequence),
            name,
        )
    )

    accepted: list[tuple[str, set[str]]] = []
    covered_seeds: set[str] = set()
    for seed in seeds:
        if seed in covered_seeds:
            continue
        nodes = expand_neighborhood(seed, graph, radius, max_nodes)
        branch_nodes = sum(
            1 for node in nodes if branch_excess(node, indegree, outdegree) > 0
        )
        if branch_nodes < min_branch_nodes:
            continue
        too_similar = False
        for _old_seed, old_nodes in accepted:
            denominator = max(1, min(len(nodes), len(old_nodes)))
            if len(nodes & old_nodes) / denominator > max_overlap_fraction:
                too_similar = True
                break
        if too_similar:
            covered_seeds.add(seed)
            continue
        accepted.append((seed, nodes))
        covered_seeds.update(nodes)
        if len(accepted) >= max_neighborhoods:
            break

    output: list[Neighborhood] = []
    for identifier, (seed, nodes) in enumerate(accepted, 1):
        coverages = [segments[node].coverage for node in nodes if segments[node].coverage > 0]
        output.append(
            Neighborhood(
                identifier=identifier,
                seed=seed,
                nodes=nodes,
                branch_nodes=sum(
                    1 for node in nodes if branch_excess(node, indegree, outdegree) > 0
                ),
                mean_coverage=sum(coverages) / len(coverages) if coverages else 0.0,
            )
        )
    return output


def build_signature_index(
    neighborhoods: list[Neighborhood],
    segments: dict[str, Segment],
    seed_k: int,
    stride: int,
) -> dict[str, int]:
    memberships: dict[str, set[int]] = defaultdict(set)
    for neighborhood in neighborhoods:
        local: set[str] = set()
        for node in neighborhood.nodes:
            local.update(canonical_kmers(segments[node].sequence, seed_k, stride))
        for kmer in local:
            memberships[kmer].add(neighborhood.identifier)
    return {
        kmer: next(iter(ids))
        for kmer, ids in memberships.items()
        if len(ids) == 1
    }


def route_pairs(
    read1: Path,
    read2: Path,
    output_dir: Path,
    neighborhoods: list[Neighborhood],
    signatures: dict[str, int],
    seed_k: int,
    stride: int,
    min_hits: int,
    hit_margin: int,
) -> int:
    handles: dict[int, tuple[TextIO, TextIO]] = {}
    by_id = {item.identifier: item for item in neighborhoods}
    for item in neighborhoods:
        directory = output_dir / f"neighborhood_{item.identifier:03d}" / "reads"
        directory.mkdir(parents=True, exist_ok=True)
        handles[item.identifier] = (
            gzip.open(directory / "R1.fastq.gz", "wt"),
            gzip.open(directory / "R2.fastq.gz", "wt"),
        )

    routed = 0
    left_iter = fastq_records(read1)
    right_iter = fastq_records(read2)
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

        scores: Counter[int] = Counter()
        for sequence in (left[1], right[1]):
            for kmer in canonical_kmers(sequence, seed_k, stride):
                identifier = signatures.get(kmer)
                if identifier is not None:
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
        routed += 1

    for out1, out2 in handles.values():
        out1.close()
        out2.close()
    return routed


def gfa_branch_count(path: Path) -> int:
    if not path.exists():
        return 0
    _segments, _graph, indegree, outdegree = parse_gfa(path)
    return sum(
        1
        for node in set(indegree) | set(outdegree)
        if branch_excess(node, indegree, outdegree) > 0
    )


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
    subprocess.run(
        [
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
        ],
        check=True,
    )
    profile = json.loads((output_dir / "run_profile.json").read_text())
    return Candidate(
        k=k,
        output_dir=output_dir,
        bases=int(profile.get("primary_bases", 0)),
        n50=int(profile.get("primary_n50", 0)),
        contigs=int(profile.get("primary_contigs", 0)),
        branches=gfa_branch_count(output_dir / "assembly.gfa"),
    )


def choose_candidate(
    candidates: list[Candidate],
    base_k: int,
    min_base_fraction: float,
    branch_fraction: float,
    n50_gain: float,
) -> tuple[Candidate, bool, str]:
    base = next((item for item in candidates if item.k == base_k), None)
    if base is None:
        raise ValueError("base-k candidate is missing")
    if base.bases == 0:
        return base, False, "base_has_no_output"

    eligible: list[Candidate] = []
    for item in candidates:
        if item.k == base_k or item.bases == 0:
            continue
        if item.bases < math.ceil(base.bases * min_base_fraction):
            continue
        simpler = item.branches <= math.floor(base.branches * branch_fraction)
        more_contiguous = item.n50 >= math.ceil(base.n50 * n50_gain)
        recovers_more = item.bases > base.bases and item.branches < base.branches
        if simpler or more_contiguous or recovers_more:
            eligible.append(item)
    if not eligible:
        return base, False, "no_higher_k_passed_gate"

    best = max(
        eligible,
        key=lambda item: (
            item.bases - base.bases,
            base.branches - item.branches,
            item.n50 - base.n50,
            -item.k,
        ),
    )
    return best, True, "higher_k_passed_gate"


def merge_outputs(
    base_fasta: Path,
    promoted_fastas: list[Path],
    output_dir: Path,
    min_contig_length: int,
) -> Path:
    scripts = Path(__file__).resolve().parent
    raw = output_dir / "adaptive_exact_union.fasta"
    final = output_dir / "adaptive_contigs.fasta"
    subprocess.run(
        [
            sys.executable,
            str(scripts / "merge_fasta_unique.py"),
            str(raw),
            str(base_fasta),
            *map(str, promoted_fastas),
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
    parser.add_argument("--routing-stride", type=int, default=2)
    parser.add_argument("--signature-stride", type=int, default=2)
    parser.add_argument("--min-seed-hits", type=int, default=2)
    parser.add_argument("--seed-hit-margin", type=int, default=1)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--max-nodes", type=int, default=64)
    parser.add_argument("--max-neighborhoods", type=int, default=16)
    parser.add_argument("--min-branch-nodes", type=int, default=2)
    parser.add_argument("--max-overlap-fraction", type=float, default=0.35)
    parser.add_argument("--min-pairs", type=int, default=30)
    parser.add_argument("--min-candidate-base-fraction", type=float, default=0.90)
    parser.add_argument("--branch-fraction", type=float, default=0.80)
    parser.add_argument("--n50-gain", type=float, default=1.15)
    parser.add_argument("--min-contig-length", type=int, default=200)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    if not args.bridgeasm.exists():
        raise FileNotFoundError(args.bridgeasm)
    if not 0.0 < args.min_candidate_base_fraction <= 1.0:
        raise ValueError("min-candidate-base-fraction must be in (0, 1]")
    if not 0.0 < args.branch_fraction <= 1.0:
        raise ValueError("branch-fraction must be in (0, 1]")
    if args.n50_gain < 1.0:
        raise ValueError("n50-gain must be >= 1")
    if not 0.0 <= args.max_overlap_fraction <= 1.0:
        raise ValueError("max-overlap-fraction must be in 0..1")

    base_gfa = args.base_dir / "assembly.gfa"
    base_fasta = args.base_dir / "primary_contigs.fasta"
    if not base_gfa.exists() or not base_fasta.exists():
        raise FileNotFoundError("base-dir must contain assembly.gfa and primary_contigs.fasta")

    candidate_ks = {args.base_k}
    candidate_ks.update(
        int(value) for value in args.candidate_k.split(",") if value.strip()
    )
    candidate_ks = {
        value for value in candidate_ks if value > 0 and value % 2 == 1
    }

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
        args.min_branch_nodes,
        args.max_overlap_fraction,
    )

    if not neighborhoods:
        final = merge_outputs(base_fasta, [], args.output, args.min_contig_length)
        summary = {
            "base_k": args.base_k,
            "candidate_k": sorted(candidate_ks),
            "neighborhoods": 0,
            "routed_pairs": 0,
            "promoted": 0,
            "signature_kmers": 0,
            "final_fasta": str(final),
        }
        (args.output / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return

    signatures = build_signature_index(
        neighborhoods, segments, args.seed_k, args.signature_stride
    )
    routed = route_pairs(
        args.read1,
        args.read2,
        args.output,
        neighborhoods,
        signatures,
        args.seed_k,
        args.routing_stride,
        args.min_seed_hits,
        args.seed_hit_margin,
    )

    promoted_fastas: list[Path] = []
    choices: list[tuple[object, ...]] = []
    for item in neighborhoods:
        directory = args.output / f"neighborhood_{item.identifier:03d}"
        if item.pair_count < args.min_pairs:
            choices.append(
                (
                    item.identifier,
                    item.seed,
                    item.pair_count,
                    args.base_k,
                    0,
                    0,
                    0,
                    False,
                    "too_few_pairs",
                )
            )
            continue

        read1 = directory / "reads" / "R1.fastq.gz"
        read2 = directory / "reads" / "R2.fastq.gz"
        candidates = [
            run_bridge(
                args.bridgeasm,
                read1,
                read2,
                directory / f"k{k}",
                k,
                args.threads,
                args.min_contig_length,
            )
            for k in sorted(candidate_ks)
        ]
        chosen, promoted, reason = choose_candidate(
            candidates,
            args.base_k,
            args.min_candidate_base_fraction,
            args.branch_fraction,
            args.n50_gain,
        )
        if promoted:
            promoted_fastas.append(chosen.output_dir / "primary_contigs.fasta")
        choices.append(
            (
                item.identifier,
                item.seed,
                item.pair_count,
                chosen.k,
                chosen.bases,
                chosen.n50,
                chosen.branches,
                promoted,
                reason,
            )
        )
        with (directory / "candidate_metrics.tsv").open("w") as handle:
            handle.write(
                "k\tprimary_bases\tprimary_n50\tprimary_contigs\tbranch_nodes\n"
            )
            for candidate in candidates:
                handle.write(
                    f"{candidate.k}\t{candidate.bases}\t{candidate.n50}\t"
                    f"{candidate.contigs}\t{candidate.branches}\n"
                )

    final = merge_outputs(
        base_fasta, promoted_fastas, args.output, args.min_contig_length
    )
    with (args.output / "neighborhoods.tsv").open("w") as handle:
        handle.write(
            "neighborhood\tseed\tnodes\tbranch_nodes\tmean_coverage\tpairs\n"
        )
        for item in neighborhoods:
            handle.write(
                f"{item.identifier}\t{item.seed}\t{len(item.nodes)}\t"
                f"{item.branch_nodes}\t{item.mean_coverage:.6f}\t{item.pair_count}\n"
            )
    with (args.output / "adaptive_k_choices.tsv").open("w") as handle:
        handle.write(
            "neighborhood\tseed\tpairs\tchosen_k\tprimary_bases\tprimary_n50\t"
            "branch_nodes\tpromoted\treason\n"
        )
        for row in choices:
            handle.write("\t".join(map(str, row)) + "\n")

    summary = {
        "base_k": args.base_k,
        "candidate_k": sorted(candidate_ks),
        "neighborhoods": len(neighborhoods),
        "routed_pairs": routed,
        "promoted": sum(1 for row in choices if row[-2]),
        "signature_kmers": len(signatures),
        "final_fasta": str(final),
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
