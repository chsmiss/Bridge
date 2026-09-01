#!/usr/bin/env python3
"""Component-local k-lifting prototype for BridgeAsm.

This extends adaptive_k_local_v2 without using a reference. The global low-k
assembly stays as the recall backbone. Reads are softly routed to branchy,
covered graph neighborhoods and low-k unitigs are converted into overlapping
pseudo-reads so local higher-k rebuilds carry forward previously supported
sequence. A higher-k patch is promoted only if it preserves local bases while
reducing branching, increasing N50, or adding sequence. Promoted patches are
merged with the base assembly and conservatively stitched by reciprocal exact
overlap.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import TextIO

import adaptive_k_local_v2 as v2

_VALID = frozenset("ACGT")


def route_pairs_soft(
    read1: Path,
    read2: Path,
    output_dir: Path,
    neighborhoods: list[v2.Neighborhood],
    signatures: dict[str, int],
    seed_k: int,
    stride: int,
    min_hits: int,
    min_top_fraction: float,
    soft_second_fraction: float,
) -> tuple[int, int, dict[int, int]]:
    handles: dict[int, tuple[TextIO, TextIO]] = {}
    soft_counts: Counter[int] = Counter()
    by_id = {item.identifier: item for item in neighborhoods}
    for item in neighborhoods:
        directory = output_dir / f"neighborhood_{item.identifier:03d}" / "reads"
        directory.mkdir(parents=True, exist_ok=True)
        handles[item.identifier] = (
            gzip.open(directory / "R1.fastq.gz", "wt"),
            gzip.open(directory / "R2.fastq.gz", "wt"),
        )

    routed_pairs = 0
    assignments = 0
    left_iter = v2.fastq_records(read1)
    right_iter = v2.fastq_records(read2)
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
        total_hits = 0
        for sequence in (left[1], right[1]):
            for kmer in v2.canonical_kmers(sequence, seed_k, stride):
                identifier = signatures.get(kmer)
                if identifier is not None:
                    scores[identifier] += 1
                    total_hits += 1
        ranked = scores.most_common(2)
        if not ranked or ranked[0][1] < min_hits or total_hits == 0:
            continue
        if ranked[0][1] / total_hits < min_top_fraction:
            continue

        selected = [ranked[0][0]]
        if (
            len(ranked) > 1
            and ranked[1][1] >= min_hits
            and ranked[1][1] >= math.ceil(ranked[0][1] * soft_second_fraction)
        ):
            selected.append(ranked[1][0])

        routed_pairs += 1
        for identifier in selected:
            out1, out2 = handles[identifier]
            v2.write_record(out1, left)
            v2.write_record(out2, right)
            by_id[identifier].pair_count += 1
            if len(selected) > 1:
                soft_counts[identifier] += 1
            assignments += 1

    for out1, out2 in handles.values():
        out1.close()
        out2.close()
    return routed_pairs, assignments, dict(soft_counts)


def make_carry_fastq(
    item: v2.Neighborhood,
    segments: dict[str, v2.Segment],
    pair_r1: Path,
    pair_r2: Path,
    output: Path,
    read_length: int,
    stride: int,
) -> int:
    """Combine real mates as single reads plus low-k unitig pseudo-reads."""
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(output, "wt") as handle:
        for mate_index, path in enumerate((pair_r1, pair_r2), 1):
            for index, record in enumerate(v2.fastq_records(path), 1):
                v2.write_record(
                    handle,
                    (f"@local_real_{mate_index}_{index}", record[1], "+", record[3]),
                )
                count += 1
        for node in sorted(item.nodes):
            sequence = segments[node].sequence
            if len(sequence) < read_length:
                continue
            starts = list(range(0, len(sequence) - read_length + 1, stride))
            last = len(sequence) - read_length
            if not starts or starts[-1] != last:
                starts.append(last)
            for start in starts:
                window = sequence[start : start + read_length]
                if set(window) - _VALID:
                    continue
                v2.write_record(
                    handle,
                    (
                        f"@carry_{item.identifier}_{node}_{start}",
                        window,
                        "+",
                        "I" * len(window),
                    ),
                )
                count += 1
    return count


def run_bridge_single(
    bridgeasm: Path,
    reads: Path,
    output_dir: Path,
    k: int,
    threads: int,
    min_contig_length: int,
) -> v2.Candidate:
    output_dir.mkdir(parents=True, exist_ok=True)
    mercy_span = 16 if k <= 31 else 8
    subprocess.run(
        [
            str(bridgeasm),
            "assemble",
            "-1",
            str(reads),
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
            "3",
            "--primary-dominance",
            "0.70",
            "--threaded-path-cover",
            "--min-contig-length",
            str(min_contig_length),
            "--threads",
            str(threads),
        ],
        check=True,
    )
    profile = json.loads((output_dir / "run_profile.json").read_text())
    return v2.Candidate(
        k=k,
        output_dir=output_dir,
        bases=int(profile.get("primary_bases", 0)),
        n50=int(profile.get("primary_n50", 0)),
        contigs=int(profile.get("primary_contigs", 0)),
        branches=v2.gfa_branch_count(output_dir / "assembly.gfa"),
    )


def local_candidate_ks(
    item: v2.Neighborhood,
    requested: list[int],
    base_k: int,
    read_length: int,
) -> list[int]:
    allowed = {base_k}
    for k in requested:
        if k <= base_k or k > read_length:
            continue
        if k <= 31:
            allowed.add(k)
        elif k <= 41 and (item.pair_count >= 40 or item.mean_coverage >= 2.0):
            allowed.add(k)
        elif k <= 55 and (item.pair_count >= 80 or item.mean_coverage >= 3.0):
            allowed.add(k)
        elif k <= 77 and (item.pair_count >= 160 or item.mean_coverage >= 4.5):
            allowed.add(k)
    return sorted(allowed)


def choose_candidate(
    candidates: list[v2.Candidate],
    base_k: int,
    min_base_fraction: float,
    branch_fraction: float,
    n50_gain: float,
    min_base_gain: float,
) -> tuple[v2.Candidate, bool, str]:
    base = next((item for item in candidates if item.k == base_k), None)
    if base is None:
        raise ValueError("base-k candidate is missing")
    if base.bases == 0:
        return base, False, "base_has_no_output"
    eligible: list[v2.Candidate] = []
    for item in candidates:
        if item.k == base_k or item.bases == 0:
            continue
        if item.bases < math.ceil(base.bases * min_base_fraction):
            continue
        simpler = item.branches <= math.floor(base.branches * branch_fraction)
        more_contiguous = item.n50 >= math.ceil(base.n50 * n50_gain)
        more_sequence = item.bases >= math.ceil(base.bases * min_base_gain)
        if simpler or more_contiguous or more_sequence:
            eligible.append(item)
    if not eligible:
        return base, False, "no_higher_k_passed_gate"
    best = max(
        eligible,
        key=lambda item: (
            item.bases - base.bases,
            base.branches - item.branches,
            item.n50 - base.n50,
            item.k,
        ),
    )
    return best, True, "higher_k_passed_gate"


def merge_and_stitch(
    base_fasta: Path,
    promoted_fastas: list[Path],
    output_dir: Path,
    min_contig_length: int,
    min_overlap: int,
) -> Path:
    scripts = Path(__file__).resolve().parent
    raw = output_dir / "local_lift_exact_union.fasta"
    contained = output_dir / "local_lift_contained.fasta"
    stitched = output_dir / "local_lift_stitched.fasta"
    final = output_dir / "local_lift_contigs.fasta"
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
            str(contained),
            "--min-length",
            str(min_contig_length),
            "--stats-json",
            str(output_dir / "containment_pre_stitch.json"),
            "--removed-tsv",
            str(output_dir / "contained_pre_stitch.tsv"),
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(scripts / "stitch_exact_overlaps.py"),
            str(stitched),
            str(contained),
            "--min-overlap",
            str(min_overlap),
            "--overlap-margin",
            "20",
            "--seed-length",
            str(min(31, min_overlap)),
            "--min-length",
            str(min_contig_length),
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(scripts / "filter_contained_fasta.py"),
            str(stitched),
            str(final),
            "--min-length",
            str(min_contig_length),
            "--stats-json",
            str(output_dir / "containment_post_stitch.json"),
            "--removed-tsv",
            str(output_dir / "contained_post_stitch.tsv"),
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
    parser.add_argument("--candidate-k", default="31,41,55,77")
    parser.add_argument("--seed-k", type=int, default=13)
    parser.add_argument("--routing-stride", type=int, default=1)
    parser.add_argument("--signature-stride", type=int, default=1)
    parser.add_argument("--min-seed-hits", type=int, default=2)
    parser.add_argument("--min-top-fraction", type=float, default=0.55)
    parser.add_argument("--soft-second-fraction", type=float, default=0.65)
    parser.add_argument("--radius", type=int, default=3)
    parser.add_argument("--max-nodes", type=int, default=96)
    parser.add_argument("--max-neighborhoods", type=int, default=24)
    parser.add_argument("--min-branch-nodes", type=int, default=1)
    parser.add_argument("--min-neighborhood-coverage", type=float, default=1.5)
    parser.add_argument("--max-overlap-fraction", type=float, default=0.40)
    parser.add_argument("--min-pairs", type=int, default=20)
    parser.add_argument("--carry-read-length", type=int, default=101)
    parser.add_argument("--carry-stride", type=int, default=20)
    parser.add_argument("--min-candidate-base-fraction", type=float, default=0.92)
    parser.add_argument("--branch-fraction", type=float, default=0.90)
    parser.add_argument("--n50-gain", type=float, default=1.08)
    parser.add_argument("--min-base-gain", type=float, default=1.01)
    parser.add_argument("--stitch-min-overlap", type=int, default=80)
    parser.add_argument("--min-contig-length", type=int, default=200)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    if not args.bridgeasm.exists():
        raise FileNotFoundError(args.bridgeasm)
    if not 0.0 < args.min_candidate_base_fraction <= 1.0:
        raise ValueError("min-candidate-base-fraction must be in (0, 1]")
    if not 0.0 < args.branch_fraction <= 1.0:
        raise ValueError("branch-fraction must be in (0, 1]")
    if args.n50_gain < 1.0 or args.min_base_gain < 1.0:
        raise ValueError("n50-gain and min-base-gain must be >= 1")

    base_gfa = args.base_dir / "assembly.gfa"
    base_fasta = args.base_dir / "primary_contigs.fasta"
    if not base_gfa.exists() or not base_fasta.exists():
        raise FileNotFoundError("base-dir must contain assembly.gfa and primary_contigs.fasta")
    requested = sorted(
        {
            int(value)
            for value in args.candidate_k.split(",")
            if value.strip() and int(value) > 0 and int(value) % 2 == 1
        }
    )

    args.output.mkdir(parents=True, exist_ok=True)
    segments, graph, indegree, outdegree = v2.parse_gfa(base_gfa)
    neighborhoods = v2.build_neighborhoods(
        segments,
        graph,
        indegree,
        outdegree,
        args.radius,
        args.max_nodes,
        args.max_neighborhoods * 2,
        args.min_branch_nodes,
        args.max_overlap_fraction,
    )
    neighborhoods = [
        item
        for item in neighborhoods
        if item.mean_coverage >= args.min_neighborhood_coverage
    ][: args.max_neighborhoods]

    signatures = v2.build_signature_index(
        neighborhoods, segments, args.seed_k, args.signature_stride
    )
    routed_pairs, pair_assignments, soft_counts = route_pairs_soft(
        args.read1,
        args.read2,
        args.output,
        neighborhoods,
        signatures,
        args.seed_k,
        args.routing_stride,
        args.min_seed_hits,
        args.min_top_fraction,
        args.soft_second_fraction,
    ) if neighborhoods else (0, 0, {})

    promoted_fastas: list[Path] = []
    choices: list[tuple[object, ...]] = []
    carry_counts: dict[int, int] = {}
    for item in neighborhoods:
        directory = args.output / f"neighborhood_{item.identifier:03d}"
        if item.pair_count < args.min_pairs:
            choices.append((item.identifier, item.seed, item.pair_count, soft_counts.get(item.identifier, 0), item.mean_coverage, args.base_k, 0, 0, 0, False, "too_few_pairs"))
            continue
        carry = directory / "carry_reads.fastq.gz"
        carry_counts[item.identifier] = make_carry_fastq(
            item,
            segments,
            directory / "reads" / "R1.fastq.gz",
            directory / "reads" / "R2.fastq.gz",
            carry,
            args.carry_read_length,
            args.carry_stride,
        )
        local_ks = local_candidate_ks(item, requested, args.base_k, args.carry_read_length)
        candidates = [
            run_bridge_single(
                args.bridgeasm,
                carry,
                directory / f"k{k}",
                k,
                args.threads,
                args.min_contig_length,
            )
            for k in local_ks
        ]
        chosen, promoted, reason = choose_candidate(
            candidates,
            args.base_k,
            args.min_candidate_base_fraction,
            args.branch_fraction,
            args.n50_gain,
            args.min_base_gain,
        )
        if promoted:
            promoted_fastas.append(chosen.output_dir / "primary_contigs.fasta")
        choices.append((item.identifier, item.seed, item.pair_count, soft_counts.get(item.identifier, 0), item.mean_coverage, chosen.k, chosen.bases, chosen.n50, chosen.branches, promoted, reason))
        with (directory / "candidate_metrics.tsv").open("w") as handle:
            handle.write("k\tprimary_bases\tprimary_n50\tprimary_contigs\tbranch_nodes\n")
            for candidate in candidates:
                handle.write(f"{candidate.k}\t{candidate.bases}\t{candidate.n50}\t{candidate.contigs}\t{candidate.branches}\n")

    final = merge_and_stitch(
        base_fasta,
        promoted_fastas,
        args.output,
        args.min_contig_length,
        args.stitch_min_overlap,
    )
    with (args.output / "neighborhoods.tsv").open("w") as handle:
        handle.write("neighborhood\tseed\tnodes\tbranch_nodes\tmean_coverage\tpairs\tsoft_assignments\tcarry_reads\n")
        for item in neighborhoods:
            handle.write(f"{item.identifier}\t{item.seed}\t{len(item.nodes)}\t{item.branch_nodes}\t{item.mean_coverage:.6f}\t{item.pair_count}\t{soft_counts.get(item.identifier, 0)}\t{carry_counts.get(item.identifier, 0)}\n")
    with (args.output / "local_k_choices.tsv").open("w") as handle:
        handle.write("neighborhood\tseed\tpairs\tsoft_assignments\tmean_coverage\tchosen_k\tprimary_bases\tprimary_n50\tbranch_nodes\tpromoted\treason\n")
        for row in choices:
            handle.write("\t".join(map(str, row)) + "\n")

    promoted_hist = Counter(int(row[5]) for row in choices if bool(row[-2]))
    summary = {
        "base_k": args.base_k,
        "candidate_k": requested,
        "neighborhoods": len(neighborhoods),
        "routed_pairs": routed_pairs,
        "pair_assignments": pair_assignments,
        "promoted": sum(promoted_hist.values()),
        "chosen_k_histogram": dict(sorted(promoted_hist.items())),
        "signature_kmers": len(signatures),
        "final_fasta": str(final),
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
