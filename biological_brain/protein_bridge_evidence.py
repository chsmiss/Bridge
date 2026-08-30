#!/usr/bin/env python3
"""Score existing assembly-graph links with protein-assembly continuity evidence.

The script never invents a nucleotide edge.  For every existing plus/plus GFA link it
reconstructs a short nucleotide context around the junction, translates all six frames,
and asks whether a stop-free peptide crossing the junction is supported by a protein
assembly such as Plass output.

The output is a TSV that can be consumed by bridgeasm-evidence-path.  Exact amino-acid
k-mers are used intentionally: they make the first implementation deterministic,
dependency-free, and conservative.  More sensitive aligners can later emit the same TSV
schema.
"""

from __future__ import annotations

import argparse
import collections
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Tuple


CODON_TABLE: Mapping[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class Link:
    source: str
    target: str
    overlap: int


@dataclass(frozen=True)
class CrossingPeptide:
    sequence: str
    boundary_aa: int
    strand: str
    frame: int


@dataclass(frozen=True)
class PeptideMatch:
    protein_id: str
    protein_start: int
    protein_end: int
    score: float
    ambiguity: float
    unique_kmers: int
    left_kmers: int
    right_kmers: int


@dataclass(frozen=True)
class Evidence:
    source: str
    target: str
    overlap: int
    protein_score: float
    unique_kmers: int
    left_kmers: int
    right_kmers: int
    ambiguity: float
    frame_consistency: float
    protein_id: str
    strand: str
    frame: int
    protein_start: int
    protein_end: int
    junction_peptide: str
    junction_boundary_aa: int
    breakpoint_class: str


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"{path}: sequence encountered before FASTA header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def parse_overlap(value: str) -> int:
    if value == "*":
        return 0
    if value.endswith("M") and value[:-1].isdigit():
        return int(value[:-1])
    raise ValueError(f"only simple M overlaps are supported, got {value!r}")


def read_gfa(path: Path) -> Tuple[Dict[str, str], List[Link], int]:
    segments: Dict[str, str] = {}
    links: List[Link] = []
    skipped_oriented = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("H\t") or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "S":
                if len(fields) < 3:
                    raise ValueError(f"{path}:{line_number}: malformed S line")
                if fields[2] == "*":
                    raise ValueError(f"{path}:{line_number}: segment sequence is required")
                segments[fields[1]] = fields[2].upper()
            elif fields[0] == "L":
                if len(fields) < 6:
                    raise ValueError(f"{path}:{line_number}: malformed L line")
                if fields[2] != "+" or fields[4] != "+":
                    skipped_oriented += 1
                    continue
                links.append(Link(fields[1], fields[3], parse_overlap(fields[5])))
    missing = [link for link in links if link.source not in segments or link.target not in segments]
    if missing:
        first = missing[0]
        raise ValueError(f"GFA link references an unknown segment: {first.source}->{first.target}")
    return segments, links, skipped_oriented


def translate(sequence: str, frame: int) -> str:
    amino_acids: List[str] = []
    stop = len(sequence) - 2
    for offset in range(frame, stop, 3):
        codon = sequence[offset : offset + 3]
        amino_acids.append(CODON_TABLE.get(codon, "X"))
    return "".join(amino_acids)


def crossing_peptides(
    junction: str,
    boundary_nt: int,
    min_side_aa: int,
    min_peptide_aa: int,
) -> List[CrossingPeptide]:
    candidates: List[CrossingPeptide] = []
    for strand, oriented, oriented_boundary in (
        ("+", junction, boundary_nt),
        ("-", reverse_complement(junction), len(junction) - boundary_nt),
    ):
        for frame in range(3):
            amino_acids = translate(oriented, frame)
            if not amino_acids:
                continue
            fragment_start = 0
            for fragment_end in range(len(amino_acids) + 1):
                at_end = fragment_end == len(amino_acids)
                at_stop = not at_end and amino_acids[fragment_end] == "*"
                if not at_end and not at_stop:
                    continue
                if fragment_end > fragment_start:
                    left_residues = 0
                    right_residues = 0
                    for aa_index in range(fragment_start, fragment_end):
                        codon_start = frame + aa_index * 3
                        codon_end = codon_start + 3
                        if codon_end <= oriented_boundary:
                            left_residues += 1
                        elif codon_start >= oriented_boundary:
                            right_residues += 1
                    peptide = amino_acids[fragment_start:fragment_end]
                    if (
                        len(peptide) >= min_peptide_aa
                        and left_residues >= min_side_aa
                        and right_residues >= min_side_aa
                    ):
                        boundary_aa = left_residues
                        candidates.append(
                            CrossingPeptide(
                                sequence=peptide,
                                boundary_aa=boundary_aa,
                                strand=strand,
                                frame=frame,
                            )
                        )
                fragment_start = fragment_end + 1
    return candidates


def build_protein_index(
    proteins: Mapping[str, str],
    aa_k: int,
    max_occurrences: int,
) -> Dict[str, List[Tuple[str, int]]]:
    raw: DefaultDict[str, List[Tuple[str, int]]] = collections.defaultdict(list)
    for protein_id, sequence in proteins.items():
        if len(sequence) < aa_k:
            continue
        for position in range(0, len(sequence) - aa_k + 1):
            kmer = sequence[position : position + aa_k]
            if "X" in kmer or "*" in kmer:
                continue
            raw[kmer].append((protein_id, position))
    return {
        kmer: locations
        for kmer, locations in raw.items()
        if 0 < len(locations) <= max_occurrences
    }


def classify_query_position(query_position: int, aa_k: int, boundary_aa: int) -> str:
    if query_position + aa_k <= boundary_aa:
        return "left"
    if query_position >= boundary_aa:
        return "right"
    return "cross"


def best_peptide_match(
    candidate: CrossingPeptide,
    index: Mapping[str, Sequence[Tuple[str, int]]],
    aa_k: int,
    min_side_kmers: int,
    min_total_kmers: int,
) -> Optional[PeptideMatch]:
    peptide = candidate.sequence
    if len(peptide) < aa_k:
        return None

    # A support cluster is keyed by protein and exact diagonal.  Exact diagonals are
    # conservative and work well for proteins assembled from the same read set.
    support: DefaultDict[Tuple[str, int], Dict[int, str]] = collections.defaultdict(dict)
    possible_kmers = len(peptide) - aa_k + 1
    for query_position in range(possible_kmers):
        kmer = peptide[query_position : query_position + aa_k]
        locations = index.get(kmer)
        if not locations:
            continue
        for protein_id, protein_position in locations:
            diagonal = protein_position - query_position
            support[(protein_id, diagonal)][query_position] = kmer

    if not support:
        return None

    ranked: List[Tuple[int, str, int, Dict[int, str]]] = []
    for (protein_id, diagonal), positions in support.items():
        ranked.append((len(positions), protein_id, diagonal, positions))
    ranked.sort(key=lambda item: (-item[0], item[1], item[2]))

    best_count, protein_id, diagonal, best_positions = ranked[0]
    second_count = ranked[1][0] if len(ranked) > 1 else 0
    left_positions = {
        position
        for position in best_positions
        if classify_query_position(position, aa_k, candidate.boundary_aa) == "left"
    }
    right_positions = {
        position
        for position in best_positions
        if classify_query_position(position, aa_k, candidate.boundary_aa) == "right"
    }
    left_kmers = len(left_positions)
    right_kmers = len(right_positions)
    ambiguity = second_count / max(1, best_count)

    coverage = best_count / max(1, possible_kmers)
    side_support = min(1.0, left_kmers / max(1, min_side_kmers)) * min(
        1.0, right_kmers / max(1, min_side_kmers)
    )
    amount_support = 1.0 - math.exp(-best_count / max(1.0, float(min_total_kmers)))
    raw_score = math.sqrt(max(0.0, coverage)) * side_support * amount_support
    score = max(0.0, min(1.0, raw_score * (1.0 - 0.5 * min(1.0, ambiguity))))

    # Return weak/one-sided matches as diagnostics.  The caller assigns a class and
    # downstream gating decides whether the edge is usable.
    return PeptideMatch(
        protein_id=protein_id,
        protein_start=max(0, diagonal),
        protein_end=max(0, diagonal + len(peptide)),
        score=score,
        ambiguity=max(0.0, min(1.0, ambiguity)),
        unique_kmers=best_count,
        left_kmers=left_kmers,
        right_kmers=right_kmers,
    )


def score_link(
    link: Link,
    segments: Mapping[str, str],
    protein_index: Mapping[str, Sequence[Tuple[str, int]]],
    junction_nt: int,
    aa_k: int,
    min_side_aa: int,
    min_peptide_aa: int,
    min_side_kmers: int,
    min_total_kmers: int,
    max_ambiguity: float,
) -> Evidence:
    source_sequence = segments[link.source]
    target_sequence = segments[link.target]
    overlap = min(link.overlap, len(target_sequence))
    left = source_sequence[-junction_nt:]
    right = target_sequence[overlap : overlap + junction_nt]
    junction = left + right
    boundary = len(left)

    candidates = crossing_peptides(
        junction=junction,
        boundary_nt=boundary,
        min_side_aa=min_side_aa,
        min_peptide_aa=min_peptide_aa,
    )
    if not candidates:
        return Evidence(
            source=link.source,
            target=link.target,
            overlap=link.overlap,
            protein_score=0.0,
            unique_kmers=0,
            left_kmers=0,
            right_kmers=0,
            ambiguity=1.0,
            frame_consistency=0.0,
            protein_id=".",
            strand=".",
            frame=-1,
            protein_start=-1,
            protein_end=-1,
            junction_peptide=".",
            junction_boundary_aa=-1,
            breakpoint_class="stop_or_short_at_junction",
        )

    diagnostics: List[Tuple[float, CrossingPeptide, PeptideMatch]] = []
    for candidate in candidates:
        match = best_peptide_match(
            candidate=candidate,
            index=protein_index,
            aa_k=aa_k,
            min_side_kmers=min_side_kmers,
            min_total_kmers=min_total_kmers,
        )
        if match is not None:
            diagnostics.append((match.score, candidate, match))

    if not diagnostics:
        representative = max(candidates, key=lambda item: len(item.sequence))
        return Evidence(
            source=link.source,
            target=link.target,
            overlap=link.overlap,
            protein_score=0.0,
            unique_kmers=0,
            left_kmers=0,
            right_kmers=0,
            ambiguity=1.0,
            frame_consistency=1.0,
            protein_id=".",
            strand=representative.strand,
            frame=representative.frame,
            protein_start=-1,
            protein_end=-1,
            junction_peptide=representative.sequence,
            junction_boundary_aa=representative.boundary_aa,
            breakpoint_class="no_protein_assembly_match",
        )

    diagnostics.sort(
        key=lambda item: (
            -item[0],
            item[2].ambiguity,
            -item[2].unique_kmers,
            item[2].protein_id,
        )
    )
    _, candidate, match = diagnostics[0]
    if match.left_kmers < min_side_kmers or match.right_kmers < min_side_kmers:
        breakpoint_class = "one_sided_protein_match"
    elif match.unique_kmers < min_total_kmers:
        breakpoint_class = "weak_protein_match"
    elif match.ambiguity > max_ambiguity:
        breakpoint_class = "ambiguous_homology"
    else:
        breakpoint_class = "same_orf_supported"

    usable_score = match.score if breakpoint_class == "same_orf_supported" else 0.0
    return Evidence(
        source=link.source,
        target=link.target,
        overlap=link.overlap,
        protein_score=usable_score,
        unique_kmers=match.unique_kmers,
        left_kmers=match.left_kmers,
        right_kmers=match.right_kmers,
        ambiguity=match.ambiguity,
        frame_consistency=1.0,
        protein_id=match.protein_id,
        strand=candidate.strand,
        frame=candidate.frame,
        protein_start=match.protein_start,
        protein_end=match.protein_end,
        junction_peptide=candidate.sequence,
        junction_boundary_aa=candidate.boundary_aa,
        breakpoint_class=breakpoint_class,
    )


def write_evidence(path: Path, rows: Iterable[Evidence]) -> None:
    header = [
        "source",
        "target",
        "overlap",
        "protein_score",
        "unique_kmers",
        "left_kmers",
        "right_kmers",
        "ambiguity",
        "frame_consistency",
        "protein_id",
        "strand",
        "frame",
        "protein_start",
        "protein_end",
        "junction_peptide",
        "junction_boundary_aa",
        "breakpoint_class",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            values = [
                row.source,
                row.target,
                str(row.overlap),
                f"{row.protein_score:.6f}",
                str(row.unique_kmers),
                str(row.left_kmers),
                str(row.right_kmers),
                f"{row.ambiguity:.6f}",
                f"{row.frame_consistency:.3f}",
                row.protein_id,
                row.strand,
                str(row.frame),
                str(row.protein_start),
                str(row.protein_end),
                row.junction_peptide,
                str(row.junction_boundary_aa),
                row.breakpoint_class,
            ]
            handle.write("\t".join(values) + "\n")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gfa", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True, help="Plass/Penguin-style protein assembly FASTA")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--junction-nt", type=int, default=450, help="nucleotides retained on each side")
    parser.add_argument("--aa-k", type=int, default=7)
    parser.add_argument("--min-side-aa", type=int, default=8)
    parser.add_argument("--min-peptide-aa", type=int, default=24)
    parser.add_argument("--min-side-kmers", type=int, default=2)
    parser.add_argument("--min-total-kmers", type=int, default=6)
    parser.add_argument("--max-kmer-occurrences", type=int, default=32)
    parser.add_argument("--max-ambiguity", type=float, default=0.45)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.junction_nt < 30:
        raise ValueError("--junction-nt must be at least 30")
    if args.aa_k < 3:
        raise ValueError("--aa-k must be at least 3")
    if args.min_side_aa < 1 or args.min_side_kmers < 1 or args.min_total_kmers < 1:
        raise ValueError("minimum support values must be positive")
    if not 0.0 <= args.max_ambiguity <= 1.0:
        raise ValueError("--max-ambiguity must be in [0, 1]")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        segments, links, skipped_oriented = read_gfa(args.gfa)
        proteins = dict(read_fasta(args.proteins))
        if not segments:
            raise ValueError("GFA contains no sequence-bearing segments")
        if not proteins:
            raise ValueError("protein FASTA contains no sequences")
        protein_index = build_protein_index(
            proteins=proteins,
            aa_k=args.aa_k,
            max_occurrences=args.max_kmer_occurrences,
        )
        rows = [
            score_link(
                link=link,
                segments=segments,
                protein_index=protein_index,
                junction_nt=args.junction_nt,
                aa_k=args.aa_k,
                min_side_aa=args.min_side_aa,
                min_peptide_aa=args.min_peptide_aa,
                min_side_kmers=args.min_side_kmers,
                min_total_kmers=args.min_total_kmers,
                max_ambiguity=args.max_ambiguity,
            )
            for link in links
        ]
        write_evidence(args.output, rows)
        supported = sum(row.breakpoint_class == "same_orf_supported" for row in rows)
        ambiguous = sum(row.breakpoint_class == "ambiguous_homology" for row in rows)
        print(
            f"protein evidence: segments={len(segments)} links={len(links)} "
            f"supported={supported} ambiguous={ambiguous} "
            f"proteins={len(proteins)} indexed_kmers={len(protein_index)} "
            f"skipped_oriented_links={skipped_oriented}",
            file=sys.stderr,
        )
        return 0
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
