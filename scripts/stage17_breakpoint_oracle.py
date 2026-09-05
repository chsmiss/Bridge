#!/usr/bin/env python3
"""Reference-aware breakpoint oracle for benchmark diagnosis only.

The assembler itself stays reference-free.  This script uses a known mock-
community reference after assembly to classify Stage10 uncovered intervals by
where their sequence evidence disappears:

raw reads -> k21 recall graph -> k31 resolve graph -> emitted Stage10 contigs.

Fractions are reported continuously.  The categorical label is only a compact
summary for deciding which mechanism to inspect next; it is not an optimization
gate and is never consumed by the assembler.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

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


@dataclass
class Gap:
    target: str
    species: str
    start: int
    end: int
    keys21: set[int]
    keys31: set[int]

    @property
    def length(self) -> int:
        return self.end - self.start


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
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


def fastq_sequences(path: Path) -> Iterator[str]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"truncated FASTQ: {path}")
            yield seq


def rolling_canonical_keys(seq: str, k: int) -> Iterator[int]:
    mask = (1 << (2 * k)) - 1
    high_shift = 2 * (k - 1)
    fwd = 0
    rev = 0
    valid = 0
    for ch in seq:
        value = BASE.get(ch)
        if value is None:
            fwd = rev = valid = 0
            continue
        fwd = ((fwd << 2) | value) & mask
        rev = (rev >> 2) | ((3 - value) << high_shift)
        valid += 1
        if valid >= k:
            yield min(fwd, rev)


def gfa_sequences(path: Path) -> Iterator[str]:
    with path.open() as handle:
        for raw in handle:
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] != "*":
                yield fields[2].upper()


def reference_catalog(reference_dir: Path) -> tuple[dict[str, str], dict[str, str]]:
    sequences: dict[str, str] = {}
    species_by_target: dict[str, str] = {}
    for path in sorted(reference_dir.iterdir()):
        if path.suffix.lower() not in {".fa", ".fna", ".fasta"}:
            continue
        species = path.stem
        for name, seq in fasta_records(path):
            sequences[name] = seq
            species_by_target[name] = species
    if not sequences:
        raise SystemExit(f"no reference FASTA records in {reference_dir}")
    return sequences, species_by_target


def paf_intervals(path: Path, *, min_block: int = 200, min_identity: float = 0.90) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            matches = int(fields[9])
            block = int(fields[10])
            if block < min_block or matches / max(1, block) < min_identity:
                continue
            target = fields[5]
            start = int(fields[7])
            end = int(fields[8])
            if end > start:
                intervals[target].append((start, end))
    return intervals


def merge_intervals(intervals: Iterable[tuple[int, int]], join_gap: int = 50) -> list[tuple[int, int]]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        cur = merged[-1]
        if start <= cur[1] + join_gap:
            cur[1] = max(cur[1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def uncovered_intervals(length: int, aligned: list[tuple[int, int]], min_gap: int = 200) -> list[tuple[int, int]]:
    if not aligned:
        return [(0, length)] if length >= min_gap else []
    gaps: list[tuple[int, int]] = []
    cursor = 0
    for start, end in aligned:
        if start - cursor >= min_gap:
            gaps.append((cursor, start))
        cursor = max(cursor, end)
    if length - cursor >= min_gap:
        gaps.append((cursor, length))
    return gaps


def sampled_sequence(seq: str, start: int, end: int, window: int = 500) -> str:
    length = end - start
    if length <= 3 * window:
        return seq[start:end]
    mid = start + length // 2
    pieces = [
        seq[start : start + window],
        seq[max(start, mid - window // 2) : min(end, mid + window // 2)],
        seq[end - window : end],
    ]
    return "N".join(pieces)


def build_gaps(
    refs: dict[str, str],
    species_by_target: dict[str, str],
    intervals: dict[str, list[tuple[int, int]]],
) -> list[Gap]:
    gaps: list[Gap] = []
    for target, seq in refs.items():
        merged = merge_intervals(intervals.get(target, []))
        for start, end in uncovered_intervals(len(seq), merged):
            sampled = sampled_sequence(seq, start, end)
            gaps.append(
                Gap(
                    target=target,
                    species=species_by_target.get(target, "unknown"),
                    start=start,
                    end=end,
                    keys21=set(rolling_canonical_keys(sampled, 21)),
                    keys31=set(rolling_canonical_keys(sampled, 31)),
                )
            )
    return gaps


def target_key_sets(gaps: list[Gap]) -> dict[int, set[int]]:
    return {
        21: set().union(*(gap.keys21 for gap in gaps)) if gaps else set(),
        31: set().union(*(gap.keys31 for gap in gaps)) if gaps else set(),
    }


def observed_target_keys(sequences: Iterable[str], k: int, targets: set[int]) -> set[int]:
    seen: set[int] = set()
    if not targets:
        return seen
    for seq in sequences:
        for key in rolling_canonical_keys(seq, k):
            if key in targets:
                seen.add(key)
        if len(seen) == len(targets):
            break
    return seen


def read_target_keys(read1: Path, read2: Path, targets: dict[int, set[int]]) -> dict[int, set[int]]:
    seen = {21: set(), 31: set()}
    for path in (read1, read2):
        for seq in fastq_sequences(path):
            for k in (21, 31):
                wanted = targets[k]
                if not wanted or len(seen[k]) == len(wanted):
                    continue
                for key in rolling_canonical_keys(seq, k):
                    if key in wanted:
                        seen[k].add(key)
    return seen


def reference_target_counts(refs: dict[str, str], targets31: set[int]) -> Counter[int]:
    counts: Counter[int] = Counter()
    for seq in refs.values():
        local = set(rolling_canonical_keys(seq, 31))
        for key in local:
            if key in targets31:
                counts[key] += 1
    return counts


def fraction(keys: set[int], observed: set[int]) -> float:
    return len(keys & observed) / len(keys) if keys else 0.0


def classify(
    raw21: float,
    raw31: float,
    graph21: float,
    graph31: float,
    emitted31: float,
) -> str:
    if raw21 < 0.15 and raw31 < 0.10:
        return "sampling_or_no_read_support"
    if graph21 < 0.20:
        return "low_k_graph_loss"
    if graph31 < 0.35 and graph31 + 0.10 < graph21:
        return "k21_to_k31_propagation_loss"
    if emitted31 < 0.35 and emitted31 + 0.10 < graph31:
        return "emission_or_path_selection_loss"
    if emitted31 < 0.60:
        return "fragmentation_or_unresolved_graph"
    return "alignment_repeat_or_postprocess"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-dir", type=Path, required=True)
    ap.add_argument("--paf", type=Path, required=True)
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--k21-gfa", type=Path, required=True)
    ap.add_argument("--k31-gfa", type=Path, required=True)
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--output-prefix", type=Path, required=True)
    args = ap.parse_args()

    refs, species_by_target = reference_catalog(args.reference_dir)
    intervals = paf_intervals(args.paf)
    gaps = build_gaps(refs, species_by_target, intervals)
    targets = target_key_sets(gaps)
    raw_seen = read_target_keys(args.read1, args.read2, targets)
    graph21_seen = observed_target_keys(gfa_sequences(args.k21_gfa), 21, targets[21])
    graph31_seen = observed_target_keys(gfa_sequences(args.k31_gfa), 31, targets[31])
    emitted31_seen = observed_target_keys(
        (seq for _name, seq in fasta_records(args.baseline)), 31, targets[31]
    )
    ref_counts31 = reference_target_counts(refs, targets[31])

    rows: list[dict[str, object]] = []
    summary: dict[tuple[str, str], dict[str, float]] = defaultdict(
        lambda: {
            "gaps": 0,
            "gap_bases": 0,
            "raw21_weighted": 0.0,
            "raw31_weighted": 0.0,
            "graph21_weighted": 0.0,
            "graph31_weighted": 0.0,
            "emitted31_weighted": 0.0,
        }
    )
    for gap in gaps:
        unique31 = {key for key in gap.keys31 if ref_counts31.get(key, 0) == 1}
        eval31 = unique31 if len(unique31) >= 10 else gap.keys31
        raw21 = fraction(gap.keys21, raw_seen[21])
        raw31 = fraction(eval31, raw_seen[31])
        graph21 = fraction(gap.keys21, graph21_seen)
        graph31 = fraction(eval31, graph31_seen)
        emitted31 = fraction(eval31, emitted31_seen)
        category = classify(raw21, raw31, graph21, graph31, emitted31)
        rows.append(
            {
                "species": gap.species,
                "target": gap.target,
                "start": gap.start,
                "end": gap.end,
                "gap_bases": gap.length,
                "sampled_k21": len(gap.keys21),
                "sampled_k31": len(gap.keys31),
                "unique_sampled_k31": len(unique31),
                "raw21_fraction": raw21,
                "raw31_fraction": raw31,
                "k21_graph_fraction": graph21,
                "k31_graph_fraction": graph31,
                "stage10_fraction": emitted31,
                "category": category,
            }
        )
        bucket = summary[(gap.species, category)]
        bucket["gaps"] += 1
        bucket["gap_bases"] += gap.length
        for name, value in (
            ("raw21_weighted", raw21),
            ("raw31_weighted", raw31),
            ("graph21_weighted", graph21),
            ("graph31_weighted", graph31),
            ("emitted31_weighted", emitted31),
        ):
            bucket[name] += value * gap.length

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    gap_path = prefix.with_name(prefix.name + "_gaps.tsv")
    summary_path = prefix.with_name(prefix.name + "_summary.tsv")
    json_path = prefix.with_name(prefix.name + "_stats.json")

    fields = [
        "species",
        "target",
        "start",
        "end",
        "gap_bases",
        "sampled_k21",
        "sampled_k31",
        "unique_sampled_k31",
        "raw21_fraction",
        "raw31_fraction",
        "k21_graph_fraction",
        "k31_graph_fraction",
        "stage10_fraction",
        "category",
    ]
    with gap_path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in (
                "raw21_fraction",
                "raw31_fraction",
                "k21_graph_fraction",
                "k31_graph_fraction",
                "stage10_fraction",
            ):
                out[key] = f"{float(out[key]):.6f}"
            writer.writerow(out)

    with summary_path.open("w") as handle:
        handle.write(
            "species\tcategory\tgaps\tgap_bases\traw21_fraction\traw31_fraction\t"
            "k21_graph_fraction\tk31_graph_fraction\tstage10_fraction\n"
        )
        for (species, category), values in sorted(
            summary.items(), key=lambda item: (item[0][0], -item[1]["gap_bases"], item[0][1])
        ):
            bases = max(1.0, values["gap_bases"])
            handle.write(
                f"{species}\t{category}\t{int(values['gaps'])}\t{int(values['gap_bases'])}\t"
                f"{values['raw21_weighted']/bases:.6f}\t"
                f"{values['raw31_weighted']/bases:.6f}\t"
                f"{values['graph21_weighted']/bases:.6f}\t"
                f"{values['graph31_weighted']/bases:.6f}\t"
                f"{values['emitted31_weighted']/bases:.6f}\n"
            )

    category_bases: Counter[str] = Counter()
    species_bases: Counter[str] = Counter()
    for row in rows:
        category_bases[str(row["category"])] += int(row["gap_bases"])
        species_bases[str(row["species"])] += int(row["gap_bases"])
    stats = {
        "gaps": len(rows),
        "gap_bases": sum(int(row["gap_bases"]) for row in rows),
        "target_k21": len(targets[21]),
        "target_k31": len(targets[31]),
        "raw_seen_k21": len(raw_seen[21]),
        "raw_seen_k31": len(raw_seen[31]),
        "graph21_seen": len(graph21_seen),
        "graph31_seen": len(graph31_seen),
        "stage10_seen": len(emitted31_seen),
        "category_gap_bases": dict(category_bases),
        "species_gap_bases": dict(species_bases),
    }
    json_path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(summary_path.read_text())
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
