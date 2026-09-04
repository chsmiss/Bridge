#!/usr/bin/env python3
"""Create dense all-vs-all Biological Brain anchor pairs only inside focus bins.

The first sparse probe intentionally uses cheap-feature-diverse anchors to detect whether a
bin is suspicious. Once a bin is promoted, using the same cheap geometry to choose more
anchors can be biased precisely when the hidden genomes have indistinguishable coverage or
composition. This second-stage sampler therefore chooses anchors by a deterministic hash of
the sequence itself. It is truth-free, independent of contig names, and approximately
random with respect to the cheap features that caused the ambiguity.

Only focus-bin anchors are emitted, so increasing from 12 to 48 or 64 anchors does not turn
DNABERT-S into a global binner.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import bridgebin_candidate_pairs as candidates


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs", type=Path, required=True)
    p.add_argument("--assignments", type=Path, required=True)
    p.add_argument("--focus-bins", type=Path, required=True)
    p.add_argument("--coverage", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--anchors-per-bin", type=int, default=48)
    p.add_argument("--min-length", type=int, default=1500)
    p.add_argument(
        "--long-core",
        type=int,
        default=4,
        help="retain this many longest contigs before deterministic sequence-hash sampling",
    )
    return p.parse_args(argv)


def sequence_hash(sequence: str) -> bytes:
    # Sequence rather than contig ID keeps the selection independent of benchmark naming.
    return hashlib.blake2b(sequence.encode("ascii", errors="ignore"), digest_size=16).digest()


def read_sequences(path: Path, min_length: int) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for contig, sequence in candidates.read_fasta(path):
        if len(sequence) >= min_length:
            out[contig] = sequence
    return out


def choose(
    members: Sequence[candidates.Feature],
    sequences: Dict[str, str],
    count: int,
    long_core: int,
) -> List[candidates.Feature]:
    if count <= 0 or not members:
        return []
    count = min(count, len(members))
    ordered = sorted(members, key=lambda feature: (-feature.length, feature.contig))
    core_count = min(max(0, long_core), count)
    selected = list(ordered[:core_count])
    selected_ids = {feature.contig for feature in selected}
    remaining = [feature for feature in members if feature.contig not in selected_ids]
    remaining.sort(
        key=lambda feature: (
            sequence_hash(sequences.get(feature.contig, "")),
            -feature.length,
            feature.contig,
        )
    )
    selected.extend(remaining[: count - len(selected)])
    return selected


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.anchors_per_bin < 2:
        raise SystemExit("--anchors-per-bin must be >= 2")
    if args.min_length < 1 or args.long_core < 0:
        raise SystemExit("length/core parameters must be non-negative")

    focus = candidates.read_focus_bins(args.focus_bins) or set()
    if not focus:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
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
        print("bridgebin-focus-anchors: focus_bins=0 anchors=0 pairs=0")
        return 0

    coverage = candidates.read_coverage(args.coverage)
    assignments = candidates.read_assignments(args.assignments)
    sequences = read_sequences(args.contigs, args.min_length)
    features = candidates.build_features(args.contigs, coverage, args.min_length)

    by_bin: Dict[str, List[candidates.Feature]] = {bin_id: [] for bin_id in focus}
    for contig, feature in features.items():
        bin_id = assignments.get(contig)
        if bin_id in by_bin:
            by_bin[bin_id].append(feature)

    records: Dict[Tuple[str, str], candidates.PairRecord] = {}
    anchor_total = 0
    per_bin: List[Tuple[str, int]] = []
    for bin_id in sorted(by_bin):
        anchors = choose(
            by_bin[bin_id], sequences, args.anchors_per_bin, args.long_core
        )
        anchor_total += len(anchors)
        per_bin.append((bin_id, len(anchors)))
        for left_index in range(len(anchors)):
            for right_index in range(left_index + 1, len(anchors)):
                candidates.add_pair(
                    records,
                    anchors[left_index],
                    anchors[right_index],
                    "focus_dense_anchor",
                )

    candidates.write_records(args.output, records, 0)
    detail = ",".join(f"{bin_id}:{count}" for bin_id, count in per_bin)
    print(
        f"bridgebin-focus-anchors: focus_bins={len(focus)} anchors={anchor_total} "
        f"pairs={len(records)} per_bin={detail}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
