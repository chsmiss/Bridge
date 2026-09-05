#!/usr/bin/env python3
"""Estimate read-level upper bounds on recoverable genome fraction from PAF alignments.

The script deliberately reports several ceilings rather than pretending that a
single mapping policy is a mathematical assembly limit:

* strict_primary: unique/high-quality primary read alignments;
* relaxed_primary: all sufficiently long, high-identity primary alignments;
* permissive_all: primary + secondary alignments meeting the relaxed sequence
  criteria (an intentionally generous, not simultaneously realizable bound);
* *_islands_ge_N: reference bases in merged covered intervals at least N bp long.

GF is computed against the total concatenated reference length, matching the
combined-reference intuition used in the Zymo MetaQUAST benchmarks.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    length = 0
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                if name is None:
                    raise ValueError("FASTA sequence before first header")
                length += len(line)
    if name is not None:
        lengths[name] = length
    if not lengths:
        raise ValueError("empty reference FASTA")
    return lengths


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[list[int]] = []
    for start, end in intervals:
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def parse_tp(optional: list[str]) -> str | None:
    for field in optional:
        if field.startswith("tp:A:"):
            return field[5:]
    return None


def summarize(
    refs: dict[str, int],
    intervals: dict[str, list[tuple[int, int]]],
    min_island: int,
) -> dict[str, object]:
    total_ref = sum(refs.values())
    total_covered = 0
    total_island = 0
    per_ref: dict[str, dict[str, float | int]] = {}
    for name, ref_len in refs.items():
        merged = merge_intervals(intervals.get(name, []))
        covered = sum(end - start for start, end in merged)
        island = sum(end - start for start, end in merged if end - start >= min_island)
        total_covered += covered
        total_island += island
        per_ref[name] = {
            "reference_bp": ref_len,
            "covered_bp": covered,
            "covered_fraction": covered / ref_len if ref_len else 0.0,
            f"islands_ge_{min_island}_bp": island,
            f"islands_ge_{min_island}_fraction": island / ref_len if ref_len else 0.0,
            "merged_intervals": len(merged),
            f"islands_ge_{min_island}_count": sum(
                1 for start, end in merged if end - start >= min_island
            ),
        }
    return {
        "total_reference_bp": total_ref,
        "covered_bp": total_covered,
        "genome_fraction_percent": 100.0 * total_covered / total_ref,
        f"islands_ge_{min_island}_bp": total_island,
        f"islands_ge_{min_island}_genome_fraction_percent": 100.0
        * total_island
        / total_ref,
        "per_reference": per_ref,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--tsv", type=Path, required=True)
    parser.add_argument("--min-island", type=int, default=200)
    parser.add_argument("--strict-identity", type=float, default=0.95)
    parser.add_argument("--strict-query-fraction", type=float, default=0.90)
    parser.add_argument("--strict-mapq", type=int, default=20)
    parser.add_argument("--relaxed-identity", type=float, default=0.90)
    parser.add_argument("--relaxed-query-fraction", type=float, default=0.80)
    args = parser.parse_args()

    refs = fasta_lengths(args.reference)
    interval_sets: dict[str, dict[str, list[tuple[int, int]]]] = {
        "strict_primary": defaultdict(list),
        "relaxed_primary": defaultdict(list),
        "permissive_all": defaultdict(list),
    }
    alignment_counts = defaultdict(int)

    with args.paf.open() as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            target = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            block_len = int(fields[10])
            mapq = int(fields[11])
            if target not in refs or qlen <= 0 or block_len <= 0:
                continue
            query_fraction = (qend - qstart) / qlen
            identity = matches / block_len
            tp = parse_tp(fields[12:])
            primary = tp != "S"

            relaxed = (
                identity >= args.relaxed_identity
                and query_fraction >= args.relaxed_query_fraction
            )
            strict = (
                primary
                and identity >= args.strict_identity
                and query_fraction >= args.strict_query_fraction
                and mapq >= args.strict_mapq
            )
            if relaxed:
                interval_sets["permissive_all"][target].append((tstart, tend))
                alignment_counts["permissive_all"] += 1
                if primary:
                    interval_sets["relaxed_primary"][target].append((tstart, tend))
                    alignment_counts["relaxed_primary"] += 1
            if strict:
                interval_sets["strict_primary"][target].append((tstart, tend))
                alignment_counts["strict_primary"] += 1

    output: dict[str, object] = {
        "parameters": {
            "min_island": args.min_island,
            "strict_identity": args.strict_identity,
            "strict_query_fraction": args.strict_query_fraction,
            "strict_mapq": args.strict_mapq,
            "relaxed_identity": args.relaxed_identity,
            "relaxed_query_fraction": args.relaxed_query_fraction,
        },
        "alignment_counts": dict(alignment_counts),
    }
    for label, intervals in interval_sets.items():
        output[label] = summarize(refs, intervals, args.min_island)

    args.json.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    with args.tsv.open("w") as handle:
        handle.write("oracle\talignments\tcovered_bp\tgf_percent\tisland_bp\tisland_gf_percent\n")
        for label in ("strict_primary", "relaxed_primary", "permissive_all"):
            summary = output[label]
            assert isinstance(summary, dict)
            handle.write(
                f"{label}\t{alignment_counts[label]}\t{summary['covered_bp']}\t"
                f"{summary['genome_fraction_percent']:.6f}\t"
                f"{summary[f'islands_ge_{args.min_island}_bp']}\t"
                f"{summary[f'islands_ge_{args.min_island}_genome_fraction_percent']:.6f}\n"
            )

    print(args.tsv.read_text(), end="")


if __name__ == "__main__":
    main()
