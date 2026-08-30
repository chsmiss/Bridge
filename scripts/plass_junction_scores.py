#!/usr/bin/env python3
"""Score candidate nucleotide joins with same-protein evidence from a PLASS catalog.

The script consumes:

* the immutable nucleotide-backbone FASTA,
* Prodigal proteins predicted from that backbone,
* MMseqs alignments from those proteins to a PLASS assembly, and
* the candidate-edge TSV emitted by ``bridgeasm-proteinguide``.

Only a positive, orientation-consistent split-protein match contributes a score.
Missing protein evidence is neutral and never creates sequence or a graph edge.
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Orf:
    query: str
    contig: str
    start: int
    end: int
    strand: int
    aa_length: int


@dataclass(frozen=True)
class Hit:
    query: str
    target: str
    identity: float
    alignment_aa: int
    target_start: int
    target_end: int
    target_length: int
    bits: float


@dataclass(frozen=True)
class OrientedOrf:
    orf: Orf
    oriented_start: int
    oriented_end: int
    oriented_strand: int
    terminal_distance: int


@dataclass(frozen=True)
class Support:
    source_orf: str
    target_orf: str
    protein: str
    sense: str
    protein_gap: int
    identity: float
    protein_coverage: float
    terminal_distance: int
    score: float


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    current = 0
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in lengths:
                        raise ValueError(f"duplicate FASTA name: {name}")
                    lengths[name] = current
                name = line[1:].split()[0]
                if not name:
                    raise ValueError("empty FASTA header")
                current = 0
            else:
                if name is None:
                    raise ValueError("sequence before first FASTA header")
                current += len(line)
    if name is not None:
        if name in lengths:
            raise ValueError(f"duplicate FASTA name: {name}")
        lengths[name] = current
    return lengths


def infer_contig(query: str, contigs: set[str]) -> str | None:
    if query in contigs:
        return query
    parts = query.split("_")
    for end in range(len(parts) - 1, 0, -1):
        candidate = "_".join(parts[:end])
        if candidate in contigs:
            return candidate
    return None


def read_prodigal_proteins(path: Path, contig_lengths: dict[str, int]) -> dict[str, Orf]:
    orfs: dict[str, Orf] = {}
    contigs = set(contig_lengths)
    query: str | None = None
    coordinates: tuple[int, int, int] | None = None
    aa_length = 0

    def commit() -> None:
        nonlocal query, coordinates, aa_length
        if query is None:
            return
        if coordinates is None:
            raise ValueError(f"Prodigal header lacks coordinates: {query}")
        contig = infer_contig(query, contigs)
        if contig is None:
            raise ValueError(f"cannot map Prodigal protein {query!r} to a backbone contig")
        start, end, strand = coordinates
        if start < 1 or end < start or end > contig_lengths[contig]:
            raise ValueError(f"invalid Prodigal coordinates for {query}: {start}-{end}")
        orfs[query] = Orf(query, contig, start, end, strand, aa_length)
        query = None
        coordinates = None
        aa_length = 0

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                fields = [field.strip() for field in line[1:].split("#")]
                query = fields[0].split()[0]
                if len(fields) < 4:
                    raise ValueError(f"unexpected Prodigal FASTA header: {line}")
                coordinates = (int(fields[1]), int(fields[2]), int(fields[3]))
                if coordinates[2] not in (-1, 1):
                    raise ValueError(f"invalid Prodigal strand for {query}: {coordinates[2]}")
            else:
                if query is None:
                    raise ValueError("protein sequence before first FASTA header")
                aa_length += len(line.replace("*", ""))
    commit()
    return orfs


def normalize_identity(value: float) -> float:
    return value / 100.0 if value > 1.0 else value


def read_hits(
    path: Path,
    orfs: dict[str, Orf],
    min_identity: float,
    min_alignment_aa: int,
    min_bits: float,
    max_hits_per_orf: int,
) -> dict[str, list[Hit]]:
    hits: dict[str, list[Hit]] = {}
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(f"alignment line {line_number} has fewer than 12 columns")
            query, target = fields[0], fields[1]
            if query not in orfs:
                continue
            identity = normalize_identity(float(fields[2]))
            alignment_aa = int(fields[3])
            target_start = min(int(fields[7]), int(fields[8])) - 1
            target_end = max(int(fields[7]), int(fields[8]))
            target_length = int(fields[9])
            bits = float(fields[11])
            if (
                identity < min_identity
                or alignment_aa < min_alignment_aa
                or bits < min_bits
                or target_start < 0
                or target_end <= target_start
                or target_end > target_length
            ):
                continue
            hits.setdefault(query, []).append(
                Hit(
                    query=query,
                    target=target,
                    identity=identity,
                    alignment_aa=alignment_aa,
                    target_start=target_start,
                    target_end=target_end,
                    target_length=target_length,
                    bits=bits,
                )
            )
    for query, query_hits in hits.items():
        query_hits.sort(key=lambda hit: (-hit.bits, -hit.identity, hit.target))
        hits[query] = query_hits[:max_hits_per_orf]
    return hits


def parse_oriented_label(label: str) -> tuple[str, bool]:
    if label.endswith("+"):
        return label[:-1], False
    if label.endswith("-"):
        return label[:-1], True
    raise ValueError(f"edge endpoint lacks orientation suffix: {label}")


def orient_orf(orf: Orf, contig_length: int, reverse: bool, source: bool) -> OrientedOrf:
    start0 = orf.start - 1
    end0 = orf.end
    if reverse:
        oriented_start = contig_length - end0
        oriented_end = contig_length - start0
        oriented_strand = -orf.strand
    else:
        oriented_start = start0
        oriented_end = end0
        oriented_strand = orf.strand
    terminal_distance = contig_length - oriented_end if source else oriented_start
    return OrientedOrf(
        orf=orf,
        oriented_start=oriented_start,
        oriented_end=oriented_end,
        oriented_strand=oriented_strand,
        terminal_distance=terminal_distance,
    )


def interval_union_length(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def support_score(
    source: OrientedOrf,
    target: OrientedOrf,
    source_hit: Hit,
    target_hit: Hit,
    sense: int,
    protein_gap: int,
    end_window: int,
    max_protein_gap: int,
    max_protein_overlap: int,
) -> Support:
    protein_coverage = interval_union_length(
        [
            (source_hit.target_start, source_hit.target_end),
            (target_hit.target_start, target_hit.target_end),
        ]
    ) / max(1, source_hit.target_length)
    identity = min(source_hit.identity, target_hit.identity)
    terminal_distance = max(source.terminal_distance, target.terminal_distance)
    terminal_score = max(0.0, 1.0 - terminal_distance / max(1, end_window))
    if protein_gap >= 0:
        continuity = max(0.0, 1.0 - protein_gap / max(1, max_protein_gap))
    else:
        continuity = max(0.0, 1.0 - (-protein_gap) / max(1, max_protein_overlap))
    aligned_fraction = min(
        1.0,
        (source_hit.alignment_aa + target_hit.alignment_aa)
        / max(1.0, source_hit.target_length * 0.5),
    )
    score = (
        0.45 * identity
        + 0.20 * continuity
        + 0.15 * terminal_score
        + 0.10 * min(1.0, protein_coverage)
        + 0.10 * aligned_fraction
    )
    return Support(
        source_orf=source.orf.query,
        target_orf=target.orf.query,
        protein=source_hit.target,
        sense="forward" if sense == 1 else "reverse",
        protein_gap=protein_gap,
        identity=identity,
        protein_coverage=protein_coverage,
        terminal_distance=terminal_distance,
        score=max(0.0, min(1.0, score)),
    )


def best_support(
    source_orfs: list[OrientedOrf],
    target_orfs: list[OrientedOrf],
    hits: dict[str, list[Hit]],
    end_window: int,
    max_protein_gap: int,
    max_protein_overlap: int,
) -> Support | None:
    best: Support | None = None
    for source in source_orfs:
        if source.terminal_distance > end_window:
            continue
        for target in target_orfs:
            if target.terminal_distance > end_window:
                continue
            if source.oriented_strand != target.oriented_strand:
                continue
            sense = source.oriented_strand
            source_by_protein = {hit.target: hit for hit in hits.get(source.orf.query, [])}
            target_by_protein = {hit.target: hit for hit in hits.get(target.orf.query, [])}
            for protein in source_by_protein.keys() & target_by_protein.keys():
                source_hit = source_by_protein[protein]
                target_hit = target_by_protein[protein]
                if source_hit.target_length != target_hit.target_length:
                    continue
                if sense == 1:
                    protein_gap = target_hit.target_start - source_hit.target_end
                else:
                    protein_gap = source_hit.target_start - target_hit.target_end
                if protein_gap > max_protein_gap or protein_gap < -max_protein_overlap:
                    continue
                support = support_score(
                    source,
                    target,
                    source_hit,
                    target_hit,
                    sense,
                    protein_gap,
                    end_window,
                    max_protein_gap,
                    max_protein_overlap,
                )
                if best is None or support.score > best.score:
                    best = support
    return best


def read_edges(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source", "target", "eligible", "projected_gap", "guide_bases"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"edge report is missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backbone", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--edge-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--details", type=Path)
    parser.add_argument("--end-window", type=int, default=900)
    parser.add_argument("--min-identity", type=float, default=0.50)
    parser.add_argument("--min-alignment-aa", type=int, default=12)
    parser.add_argument("--min-bits", type=float, default=30.0)
    parser.add_argument("--max-hits-per-orf", type=int, default=20)
    parser.add_argument("--max-protein-gap", type=int, default=80)
    parser.add_argument("--max-protein-overlap", type=int, default=25)
    parser.add_argument("--min-output-score", type=float, default=0.0)
    args = parser.parse_args()

    if args.end_window <= 0:
        parser.error("--end-window must be positive")
    if not 0.0 <= args.min_identity <= 1.0:
        parser.error("--min-identity must be in [0, 1]")
    if not 0.0 <= args.min_output_score <= 1.0:
        parser.error("--min-output-score must be in [0, 1]")

    contig_lengths = fasta_lengths(args.backbone)
    orfs = read_prodigal_proteins(args.proteins, contig_lengths)
    by_contig: dict[str, list[Orf]] = {}
    for orf in orfs.values():
        by_contig.setdefault(orf.contig, []).append(orf)
    hits = read_hits(
        args.alignments,
        orfs,
        args.min_identity,
        args.min_alignment_aa,
        args.min_bits,
        args.max_hits_per_orf,
    )
    edges = read_edges(args.edge_report)

    details_rows: list[dict[str, object]] = []
    output_rows: list[tuple[str, str, float]] = []
    eligible_edges = 0
    for edge in edges:
        if edge["eligible"].lower() != "true":
            continue
        eligible_edges += 1
        source_name, source_reverse = parse_oriented_label(edge["source"])
        target_name, target_reverse = parse_oriented_label(edge["target"])
        if source_name not in contig_lengths or target_name not in contig_lengths:
            raise ValueError(f"edge references unknown contig: {edge['source']} -> {edge['target']}")
        source_orfs = [
            orient_orf(orf, contig_lengths[source_name], source_reverse, source=True)
            for orf in by_contig.get(source_name, [])
        ]
        target_orfs = [
            orient_orf(orf, contig_lengths[target_name], target_reverse, source=False)
            for orf in by_contig.get(target_name, [])
        ]
        support = best_support(
            source_orfs,
            target_orfs,
            hits,
            args.end_window,
            args.max_protein_gap,
            args.max_protein_overlap,
        )
        if support is None:
            continue
        details_rows.append(
            {
                "source": edge["source"],
                "target": edge["target"],
                "score": f"{support.score:.6f}",
                "protein": support.protein,
                "source_orf": support.source_orf,
                "target_orf": support.target_orf,
                "sense": support.sense,
                "protein_gap": support.protein_gap,
                "identity": f"{support.identity:.6f}",
                "protein_coverage": f"{support.protein_coverage:.6f}",
                "terminal_distance": support.terminal_distance,
                "projected_gap": edge["projected_gap"],
                "guide_bases": edge["guide_bases"],
            }
        )
        if support.score >= args.min_output_score:
            output_rows.append((edge["source"], edge["target"], support.score))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write("source\ttarget\tscore\tdecision\tscorer\n")
        for source, target, score in sorted(output_rows):
            handle.write(f"{source}\t{target}\t{score:.6f}\tneutral\tplass_same_protein\n")

    if args.details:
        args.details.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "source",
            "target",
            "score",
            "protein",
            "source_orf",
            "target_orf",
            "sense",
            "protein_gap",
            "identity",
            "protein_coverage",
            "terminal_distance",
            "projected_gap",
            "guide_bases",
        ]
        with args.details.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(details_rows)

    print(
        "\t".join(
            [
                f"eligible_edges={eligible_edges}",
                f"supported_edges={len(details_rows)}",
                f"output_scores={len(output_rows)}",
                f"orfs={len(orfs)}",
                f"aligned_orfs={len(hits)}",
            ]
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
