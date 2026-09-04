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
    parser.add_argument(
        "--anchor-strategy",
        choices=("longest", "diverse"),
        default="longest",
        help=(
            "anchor selection policy: longest preserves the legacy behavior; diverse keeps "
            "two long cores then uses cheap-feature farthest-point sampling"
        ),
    )
    parser.add_argument("--within-neighbors", type=int, default=3)
    parser.add_argument("--within-contrast", type=int, default=2)
    parser.add_argument("--merge-bin-neighbors", type=int, default=4)
    parser.add_argument("--cross-anchor-pairs", type=int, default=8)
    parser.add_argument("--rescue-bin-neighbors", type=int, default=4)
    parser.add_argument("--rescue-anchor-pairs", type=int, default=4)
    parser.add_argument("--max-pairs", type=int, default=250000)
    parser.add_argument(
        "--focus-bins",
        type=Path,
        help=(
            "optional TSV/text file of bin IDs; per-member within-bin expansion is limited "
            "to these bins while anchor probes remain global"
        ),
    )
    parser.add_argument(
        "--skip-residual-rescue",
        action="store_true",
        help="omit residual-to-bin candidate expansion in this pass",
    )
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


def read_focus_bins(path: Optional[Path]) -> Optional[Set[str]]:
    if path is None:
        return None
    bins: Set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        first_line = handle.readline()
        if not first_line:
            return bins
        fields = first_line.rstrip("\n").split("\t")
        if fields and fields[0].strip().lower() in {"bin", "bin_id", "cluster"}:
            for raw in handle:
                value = raw.rstrip("\n").split("\t", 1)[0].strip()
                if value:
                    bins.add(value)
        else:
            value = fields[0].strip() if fields else ""
            if value:
                bins.add(value)
            for raw in handle:
                value = raw.rstrip("\n").split("\t", 1)[0].strip()
                if value:
                    bins.add(value)
    return bins


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


def cosine(left: Sequence[float], right: Sequence[float]) -> float:
    if not left or len(left) != len(right):
        return 0.0
    dot = sum(a * b for a, b in zip(left, right))
    lnorm = math.sqrt(sum(value * value for value in left))
    rnorm = math.sqrt(sum(value * value for value in right))
    if lnorm <= 1e-12 or rnorm <= 1e-12:
        return 0.0
    return max(-1.0, min(1.0, dot / (lnorm * rnorm)))


def log_coverage_similarity(left: Sequence[float], right: Sequence[float]) -> Optional[float]:
    if not left or len(left) != len(right):
        return None
    distance = sum(abs(math.log((a + 0.5) / (b + 0.5))) for a, b in zip(left, right)) / len(left)
    return math.exp(-distance / 0.85)


def cheap_similarity(left: Feature, right: Feature) -> float:
    comp = max(0.0, cosine(left.composition, right.composition))
    gc = math.exp(-abs(left.gc - right.gc) / 0.08)
    cov = log_coverage_similarity(left.coverage, right.coverage)
    weighted = 0.45 * comp + 0.05 * gc
    total = 0.50
    if cov is not None:
        weighted += 0.50 * cov
        total += 0.50
    return weighted / total


def build_features(
    contigs: Path, coverage: Dict[str, List[float]], min_length: int
) -> Dict[str, Feature]:
    result: Dict[str, Feature] = {}
    for contig, sequence in read_fasta(contigs):
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


def farthest_anchor_selection(candidates: Sequence[Feature], limit: int) -> List[Feature]:
    if limit <= 0 or not candidates:
        return []
    ordered = sorted(candidates, key=lambda item: (-item.length, item.contig))
    selected = ordered[: min(2, limit)]
    selected_names = {item.contig for item in selected}
    while len(selected) < min(limit, len(ordered)):
        remaining = [item for item in ordered if item.contig not in selected_names]
        best = min(
            remaining,
            key=lambda item: (
                max(cheap_similarity(item, anchor) for anchor in selected),
                -item.length,
                item.contig,
            ),
        )
        selected.append(best)
        selected_names.add(best.contig)
    return selected


def nearest(items: Sequence[Tuple[float, str]], limit: int, reverse: bool = True) -> List[str]:
    if limit <= 0:
        return []
    ordered = sorted(items, key=lambda item: (item[0], item[1]), reverse=reverse)
    return [name for _, name in ordered[:limit]]


def add_record(
    records: Dict[Tuple[str, str], PairRecord],
    left: Feature,
    right: Feature,
    candidate_class: str,
) -> None:
    a, b = sorted((left.contig, right.contig))
    if a == b:
        return
    key = (a, b)
    record = records.get(key)
    if record is None:
        record = PairRecord(
            left=a,
            right=b,
            coverage_similarity=log_coverage_similarity(left.coverage, right.coverage),
            composition_similarity=max(0.0, cosine(left.composition, right.composition)),
            gc_similarity=math.exp(-abs(left.gc - right.gc) / 0.08),
            classes=set(),
        )
        records[key] = record
    record.classes.add(candidate_class)


def choose_anchors(
    features: Dict[str, Feature],
    members: Sequence[str],
    limit: int,
    strategy: str,
) -> List[Feature]:
    candidates = [features[name] for name in members if name in features]
    if strategy == "diverse":
        return farthest_anchor_selection(candidates, limit)
    return sorted(candidates, key=lambda item: (-item.length, item.contig))[:limit]


def mine(
    features: Dict[str, Feature], assignments: Dict[str, Optional[str]], args: argparse.Namespace
) -> List[PairRecord]:
    bins: Dict[str, List[str]] = defaultdict(list)
    residuals: List[str] = []
    for contig, bin_name in assignments.items():
        if contig not in features:
            continue
        if bin_name is None:
            residuals.append(contig)
        else:
            bins[bin_name].append(contig)

    anchors: Dict[str, List[Feature]] = {
        bin_name: choose_anchors(features, members, args.anchors_per_bin, args.anchor_strategy)
        for bin_name, members in bins.items()
    }
    records: Dict[Tuple[str, str], PairRecord] = {}

    for bin_name, members in bins.items():
        current_anchors = anchors[bin_name]
        if len(current_anchors) < 2:
            continue
        for left_index, left in enumerate(current_anchors):
            for right in current_anchors[left_index + 1 :]:
                add_record(records, left, right, "within_bin_anchor")
        if args.probe_only:
            continue
        if args.focus_bins is not None and bin_name not in args.focus_bins:
            continue
        anchor_names = {anchor.contig for anchor in current_anchors}
        for member_name in members:
            if member_name in anchor_names:
                continue
            member = features[member_name]
            similarities = [
                (cheap_similarity(member, anchor), anchor.contig) for anchor in current_anchors
            ]
            for anchor_name in nearest(similarities, args.within_neighbors, reverse=True):
                add_record(records, member, features[anchor_name], "within_bin_member")
            for anchor_name in nearest(similarities, args.within_contrast, reverse=False):
                add_record(records, member, features[anchor_name], "within_bin_contrast")

    bin_names = sorted(anchors)
    if args.merge_bin_neighbors > 0 and args.cross_anchor_pairs > 0:
        bin_centroids: Dict[str, Feature] = {}
        for bin_name, current_anchors in anchors.items():
            if not current_anchors:
                continue
            total_length = sum(anchor.length for anchor in current_anchors)
            width = max((len(anchor.coverage) for anchor in current_anchors), default=0)
            coverage_values = [0.0] * width
            composition_values = [0.0] * 256
            gc = 0.0
            for anchor in current_anchors:
                weight = anchor.length / max(1, total_length)
                gc += anchor.gc * weight
                for index, value in enumerate(anchor.composition):
                    composition_values[index] += value * weight
                if len(anchor.coverage) == width:
                    for index, value in enumerate(anchor.coverage):
                        coverage_values[index] += value * weight
            bin_centroids[bin_name] = Feature(
                contig=bin_name,
                length=total_length,
                gc=gc,
                composition=composition_values,
                coverage=coverage_values,
            )

        for bin_name in bin_names:
            if bin_name not in bin_centroids:
                continue
            current = bin_centroids[bin_name]
            others = [
                (cheap_similarity(current, bin_centroids[other]), other)
                for other in bin_names
                if other != bin_name and other in bin_centroids
            ]
            for other in nearest(others, args.merge_bin_neighbors, reverse=True):
                left_anchors = anchors.get(bin_name, [])
                right_anchors = anchors.get(other, [])
                cross = [
                    (cheap_similarity(left, right), left, right)
                    for left in left_anchors
                    for right in right_anchors
                ]
                cross.sort(key=lambda item: (-item[0], item[1].contig, item[2].contig))
                for _, left, right in cross[: args.cross_anchor_pairs]:
                    add_record(records, left, right, "cross_bin_merge")

    if not args.probe_only and not args.skip_residual_rescue and args.rescue_bin_neighbors > 0 and args.rescue_anchor_pairs > 0:
        centroids = {
            bin_name: max(current_anchors, key=lambda item: item.length)
            for bin_name, current_anchors in anchors.items()
            if current_anchors
        }
        for residual_name in residuals:
            residual = features[residual_name]
            candidates = [
                (cheap_similarity(residual, centroid), bin_name)
                for bin_name, centroid in centroids.items()
            ]
            for bin_name in nearest(candidates, args.rescue_bin_neighbors, reverse=True):
                scores = [
                    (cheap_similarity(residual, anchor), anchor.contig)
                    for anchor in anchors[bin_name]
                ]
                for anchor_name in nearest(scores, args.rescue_anchor_pairs, reverse=True):
                    add_record(records, residual, features[anchor_name], "residual_rescue")
    return list(records.values())


def write_records(path: Path, records: Sequence[PairRecord], max_pairs: int) -> None:
    ranked = sorted(
        records,
        key=lambda record: (
            -len(record.classes),
            record.left,
            record.right,
        ),
    )[:max_pairs]
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
    if args.anchors_per_bin < 1 or args.min_length < 1:
        raise SystemExit("--anchors-per-bin and --min-length must be positive")
    for name in (
        "within_neighbors",
        "merge_bin_neighbors",
        "cross_anchor_pairs",
        "rescue_bin_neighbors",
        "rescue_anchor_pairs",
    ):
        if getattr(args, name) < 0:
            raise SystemExit(f"--{name.replace('_', '-')} must be >=0")
    if args.within_contrast < 0:
        raise SystemExit("--within-contrast must be >=0")
    coverage = read_coverage(args.coverage)
    assignments = read_assignments(args.assignments)
    args.focus_bins = read_focus_bins(args.focus_bins)
    features = build_features(args.contigs, coverage, args.min_length)
    records = mine(features, assignments, args)
    write_records(args.output, records, args.max_pairs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
