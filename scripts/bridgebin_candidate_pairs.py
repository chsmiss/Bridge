#!/usr/bin/env python3
"""Mine sparse, truth-free candidate pairs for BridgeBin Biological Brain.

The expensive same-genome model must not run on all O(N^2) contig pairs. This miner uses
only cheap evidence (coverage, canonical 4-mer composition, GC, current v2 assignments)
to achieve high *recall* of three decision classes:

1. ``within_bin_anchor`` / ``within_bin_contrast``
   probe a current bin for hidden mixtures; long anchors are compared both to their most
   compatible and deliberately least-compatible anchors so contamination can generate a
   learned hard negative.
2. ``cross_bin_merge``
   compare anchors from nearby v2 bins that may be fragments of one genome.
3. ``residual_rescue``
   compare an unbinned contig to anchors from its nearest candidate bins.

Cheap similarity only decides which pairs are sent to the Biological Brain. It never
forces a merge. Output is directly consumable by ``bridgebin_pair_head.py score``.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


@dataclass
class Feature:
    contig: str
    length: int
    gc: float
    composition: List[float]
    coverage: List[float]


@dataclass
class PairRecord:
    left: str
    right: str
    coverage_similarity: Optional[float]
    composition_similarity: float
    gc_similarity: float
    classes: Set[str]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--coverage", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-length", type=int, default=1000)
    parser.add_argument("--anchors-per-bin", type=int, default=10)
    parser.add_argument("--within-neighbors", type=int, default=3)
    parser.add_argument("--within-contrast", type=int, default=2)
    parser.add_argument("--merge-bin-neighbors", type=int, default=4)
    parser.add_argument("--cross-anchor-pairs", type=int, default=8)
    parser.add_argument("--rescue-bin-neighbors", type=int, default=4)
    parser.add_argument("--rescue-anchor-pairs", type=int, default=4)
    parser.add_argument("--max-pairs", type=int, default=250000)
    parser.add_argument(
        "--probe-only",
        action="store_true",
        help=(
            "emit only anchor-level within-bin and cross-bin pairs; skip member expansion "
            "and residual rescue for a cheap first-pass conflict probe"
        ),
    )
    return parser.parse_args(argv)


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
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
                if name is None:
                    raise ValueError(f"{path}: sequence before first FASTA header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def read_coverage(path: Optional[Path]) -> Dict[str, List[float]]:
    if path is None:
        return {}
    result: Dict[str, List[float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError(f"{path}: coverage table needs contig plus sample columns")
        width = len(header) - 1
        for line_number, fields in enumerate(reader, start=2):
            if not fields or not fields[0].strip():
                continue
            if len(fields) != width + 1:
                raise ValueError(f"{path}:{line_number}: inconsistent coverage width")
            values = [float(value) for value in fields[1:]]
            if not all(math.isfinite(value) and value >= 0.0 for value in values):
                raise ValueError(f"{path}:{line_number}: invalid coverage value")
            result[fields[0].strip()] = values
    return result


def read_assignments(path: Path) -> Dict[str, Optional[str]]:
    result: Dict[str, Optional[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "contig" not in reader.fieldnames or "bin" not in reader.fieldnames:
            raise ValueError(f"{path}: assignments must contain contig and bin columns")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            raw_bin = (row.get("bin") or "").strip()
            if not contig:
                continue
            result[contig] = None if raw_bin in {"", ".", "unbinned", "NA"} else raw_bin
    return result


def base_code(base: str) -> Optional[int]:
    return {"A": 0, "C": 1, "G": 2, "T": 3}.get(base)


def reverse_complement_code(kmer: str) -> Optional[int]:
    code = 0
    for base in reversed(kmer):
        value = base_code(base)
        if value is None:
            return None
        code = (code << 2) | (3 - value)
    return code


def forward_code(kmer: str) -> Optional[int]:
    code = 0
    for base in kmer:
        value = base_code(base)
        if value is None:
            return None
        code = (code << 2) | value
    return code


def canonical_fourmer(sequence: str) -> List[float]:
    counts = [0.0] * 256
    total = 0.0
    for start in range(0, max(0, len(sequence) - 3)):
        kmer = sequence[start : start + 4]
        forward = forward_code(kmer)
        reverse = reverse_complement_code(kmer)
        if forward is None or reverse is None:
            continue
        counts[min(forward, reverse)] += 1.0
        total += 1.0
    if total > 0.0:
        counts = [value / total for value in counts]
    return counts


def gc_fraction(sequence: str) -> float:
    valid = [base for base in sequence if base in "ACGT"]
    if not valid:
        return 0.0
    return sum(base in "GC" for base in valid) / len(valid)


def build_features(
    contig_path: Path, coverage: Dict[str, List[float]], min_length: int
) -> Dict[str, Feature]:
    result: Dict[str, Feature] = {}
    for contig, sequence in read_fasta(contig_path):
        if len(sequence) < min_length:
            continue
        result[contig] = Feature(
            contig=contig,
            length=len(sequence),
            gc=gc_fraction(sequence),
            composition=canonical_fourmer(sequence),
            coverage=coverage.get(contig, []),
        )
    return result


def hellinger(left: Sequence[float], right: Sequence[float]) -> float:
    return math.sqrt(
        0.5
        * sum((math.sqrt(a) - math.sqrt(b)) ** 2 for a, b in zip(left, right))
    )


def component_sim(left: Feature, right: Feature) -> Tuple[Optional[float], float, float]:
    composition = math.exp(-hellinger(left.composition, right.composition) / 0.30)
    gc = math.exp(-abs(left.gc - right.gc) / 0.08)
    coverage: Optional[float] = None
    if left.coverage and len(left.coverage) == len(right.coverage):
        distance = sum(
            abs(math.log((a + 0.5) / (b + 0.5)))
            for a, b in zip(left.coverage, right.coverage)
        ) / len(left.coverage)
        coverage = math.exp(-distance / 0.85)
    return coverage, composition, gc


def cheap_similarity(left: Feature, right: Feature) -> float:
    coverage, composition, gc = component_sim(left, right)
    values = [composition, gc]
    weights = [0.45, 0.15]
    if coverage is not None:
        values.append(coverage)
        weights.append(0.40)
    return sum(value * weight for value, weight in zip(values, weights)) / sum(weights)


def centroid(members: Sequence[Feature], name: str) -> Feature:
    total_bp = sum(member.length for member in members)
    if total_bp <= 0:
        raise ValueError("cannot build centroid of empty bin")
    gc = sum(member.gc * member.length for member in members) / total_bp
    composition = [0.0] * 256
    for member in members:
        weight = member.length / total_bp
        for index, value in enumerate(member.composition):
            composition[index] += value * weight
    coverage: List[float] = []
    widths = {len(member.coverage) for member in members if member.coverage}
    if len(widths) == 1:
        width = next(iter(widths))
        coverage = [0.0] * width
        covered_bp = sum(member.length for member in members if len(member.coverage) == width)
        if covered_bp > 0:
            for member in members:
                if len(member.coverage) != width:
                    continue
                weight = member.length / covered_bp
                for index, value in enumerate(member.coverage):
                    coverage[index] += value * weight
    return Feature(name, total_bp, gc, composition, coverage)


def ordered_pair(left: str, right: str) -> Tuple[str, str]:
    return (left, right) if left <= right else (right, left)


def add_pair(
    records: Dict[Tuple[str, str], PairRecord],
    left: Feature,
    right: Feature,
    candidate_class: str,
) -> None:
    if left.contig == right.contig:
        return
    a, b = (left, right) if left.contig <= right.contig else (right, left)
    key = (a.contig, b.contig)
    coverage, composition, gc = component_sim(a, b)
    existing = records.get(key)
    if existing is None:
        records[key] = PairRecord(
            left=a.contig,
            right=b.contig,
            coverage_similarity=coverage,
            composition_similarity=composition,
            gc_similarity=gc,
            classes={candidate_class},
        )
    else:
        existing.classes.add(candidate_class)


def ranked_anchors(member: Feature, anchors: Sequence[Feature]) -> List[Tuple[float, Feature]]:
    ranked = [
        (cheap_similarity(member, anchor), anchor)
        for anchor in anchors
        if anchor.contig != member.contig
    ]
    ranked.sort(key=lambda item: (-item[0], -item[1].length, item[1].contig))
    return ranked


def mine(
    features: Dict[str, Feature],
    assignments: Dict[str, Optional[str]],
    args: argparse.Namespace,
) -> Dict[Tuple[str, str], PairRecord]:
    bins: Dict[str, List[Feature]] = defaultdict(list)
    residuals: List[Feature] = []
    for contig, feature in features.items():
        bin_name = assignments.get(contig)
        if bin_name is None:
            residuals.append(feature)
        else:
            bins[bin_name].append(feature)
    for members in bins.values():
        members.sort(key=lambda feature: (-feature.length, feature.contig))

    anchors: Dict[str, List[Feature]] = {
        bin_name: members[: args.anchors_per_bin]
        for bin_name, members in bins.items()
        if members
    }
    records: Dict[Tuple[str, str], PairRecord] = {}

    # Probe every current bin for hidden mixtures. Anchor all-vs-all ensures a mixed bin
    # with two long genome components presents cross-component examples to the model.
    for bin_name, members in bins.items():
        bin_anchors = anchors.get(bin_name, [])
        for left_index in range(len(bin_anchors)):
            for right_index in range(left_index + 1, len(bin_anchors)):
                add_pair(records, bin_anchors[left_index], bin_anchors[right_index], "within_bin_anchor")
        if not args.probe_only:
            for member in members:
                ranked = ranked_anchors(member, bin_anchors)
                for _score, anchor in ranked[: args.within_neighbors]:
                    add_pair(records, member, anchor, "within_bin_neighbor")
                if args.within_contrast > 0:
                    for _score, anchor in ranked[-args.within_contrast :]:
                        add_pair(records, member, anchor, "within_bin_contrast")

    # Current bins that are close under cheap features are candidates for conservative
    # biological re-merge. Select neighbors by bin centroid, then only score anchor pairs.
    centroids: Dict[str, Feature] = {
        name: centroid(members, f"__bin__{name}") for name, members in bins.items() if members
    }
    processed_bin_pairs: Set[Tuple[str, str]] = set()
    for bin_name, center in centroids.items():
        neighbors = [
            (cheap_similarity(center, other), other_name)
            for other_name, other in centroids.items()
            if other_name != bin_name
        ]
        neighbors.sort(key=lambda item: (-item[0], item[1]))
        for _score, other_name in neighbors[: args.merge_bin_neighbors]:
            bin_pair = ordered_pair(bin_name, other_name)
            if bin_pair in processed_bin_pairs:
                continue
            processed_bin_pairs.add(bin_pair)
            possibilities = [
                (cheap_similarity(left, right), left, right)
                for left in anchors.get(bin_name, [])
                for right in anchors.get(other_name, [])
            ]
            possibilities.sort(
                key=lambda item: (-item[0], -item[1].length - item[2].length, item[1].contig, item[2].contig)
            )
            for _pair_score, left, right in possibilities[: args.cross_anchor_pairs]:
                add_pair(records, left, right, "cross_bin_merge")

    # Residual rescue is deferred until after the anchor conflict probe.  This keeps the
    # first-pass foundation-model endpoint set proportional to bin anchors, not all contigs.
    if not args.probe_only:
        for residual in sorted(residuals, key=lambda feature: (-feature.length, feature.contig)):
            neighbor_bins = [
                (cheap_similarity(residual, center), bin_name)
                for bin_name, center in centroids.items()
            ]
            neighbor_bins.sort(key=lambda item: (-item[0], item[1]))
            for _bin_score, bin_name in neighbor_bins[: args.rescue_bin_neighbors]:
                ranked = ranked_anchors(residual, anchors.get(bin_name, []))
                for _score, anchor in ranked[: args.rescue_anchor_pairs]:
                    add_pair(records, residual, anchor, "residual_rescue")

    return records


def write_records(path: Path, records: Dict[Tuple[str, str], PairRecord], max_pairs: int) -> None:
    ranked = sorted(
        records.values(),
        key=lambda record: (
            -max(
                record.composition_similarity,
                record.coverage_similarity if record.coverage_similarity is not None else 0.0,
            ),
            record.left,
            record.right,
        ),
    )
    if max_pairs > 0:
        ranked = ranked[:max_pairs]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "left",
                "right",
                "coverage_similarity",
                "composition_similarity",
                "gc_similarity",
                "physical_support",
                "candidate_class",
            ]
        )
        for record in ranked:
            writer.writerow(
                [
                    record.left,
                    record.right,
                    "."
                    if record.coverage_similarity is None
                    else f"{record.coverage_similarity:.8f}",
                    f"{record.composition_similarity:.8f}",
                    f"{record.gc_similarity:.8f}",
                    ".",
                    ";".join(sorted(record.classes)),
                ]
            )
    class_counts: Dict[str, int] = defaultdict(int)
    for record in ranked:
        for candidate_class in record.classes:
            class_counts[candidate_class] += 1
    summary = " ".join(f"{key}={value}" for key, value in sorted(class_counts.items()))
    print(f"bridgebin-candidates: pairs={len(ranked)} {summary}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    for name in (
        "anchors_per_bin",
        "within_neighbors",
        "merge_bin_neighbors",
        "cross_anchor_pairs",
        "rescue_bin_neighbors",
        "rescue_anchor_pairs",
    ):
        if getattr(args, name) < 1:
            raise SystemExit(f"--{name.replace('_', '-')} must be positive")
    if args.within_contrast < 0 or args.min_length < 1:
        raise SystemExit("--within-contrast must be >=0 and --min-length must be positive")
    coverage = read_coverage(args.coverage)
    assignments = read_assignments(args.assignments)
    features = build_features(args.contigs, coverage, args.min_length)
    records = mine(features, assignments, args)
    write_records(args.output, records, args.max_pairs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
