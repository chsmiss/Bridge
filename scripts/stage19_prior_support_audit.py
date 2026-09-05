#!/usr/bin/env python3
"""Audit raw physical-fragment support for validated short prior contigs.

This is reference-free diagnostic code.  Given Stage10 cross-k validated rare
seed contigs, it asks what happens to their k31/k41/k55 sequence in the raw
paired library before any synthetic projection:

* zero fragments: no exact evidence exists at that target k;
* one fragment: real sequence exists but production min-count=2 rejects it;
* two or more fragments: it should already be solid, so a later graph/path
  stage is responsible if it disappears.

The report also measures the longest consecutive runs with >=1 and >=2 physical
fragment support.  A trusted-prior singleton channel is justified only when a
substantial amount of validated sequence is in the exactly-one-fragment class;
it must never invent zero-support target-k sequence.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

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
class SeedKmers:
    seed_id: int
    name: str
    length: int
    k: int
    ordered: list[int | None]
    unique: set[int]


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


def fastq_records(path: Path) -> Iterator[tuple[str, str]]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            if not qual:
                raise ValueError(f"truncated FASTQ: {path}")
            yield seq, qual


def canonical_key(seq: str) -> int | None:
    fwd = 0
    rev = 0
    k = len(seq)
    for i, ch in enumerate(seq):
        value = BASE.get(ch)
        if value is None:
            return None
        fwd = (fwd << 2) | value
        rev |= (3 - value) << (2 * i)
    return min(fwd, rev)


def ordered_keys(seq: str, k: int) -> list[int | None]:
    if len(seq) < k:
        return []
    return [canonical_key(seq[pos : pos + k]) for pos in range(len(seq) - k + 1)]


def rolling_target_keys(seq: str, k: int, targets: set[int]) -> set[int]:
    if not targets or len(seq) < k:
        return set()
    mask = (1 << (2 * k)) - 1
    high_shift = 2 * (k - 1)
    fwd = 0
    rev = 0
    valid = 0
    found: set[int] = set()
    for ch in seq:
        value = BASE.get(ch)
        if value is None:
            fwd = rev = valid = 0
            continue
        fwd = ((fwd << 2) | value) & mask
        rev = (rev >> 2) | ((3 - value) << high_shift)
        valid += 1
        if valid >= k:
            key = min(fwd, rev)
            if key in targets:
                found.add(key)
    return found


def build_seed_kmers(seed_path: Path, ks: list[int]) -> tuple[list[SeedKmers], dict[int, set[int]], dict[int, Counter[int]]]:
    records = list(fasta_records(seed_path))
    seed_kmers: list[SeedKmers] = []
    targets: dict[int, set[int]] = {k: set() for k in ks}
    membership: dict[int, Counter[int]] = {k: Counter() for k in ks}
    for sid, (name, seq) in enumerate(records):
        for k in ks:
            ordered = ordered_keys(seq, k)
            unique = {key for key in ordered if key is not None}
            seed_kmers.append(SeedKmers(sid, name, len(seq), k, ordered, unique))
            targets[k].update(unique)
            membership[k].update(unique)
    return seed_kmers, targets, membership


def scan_fragment_support(
    read1: Path,
    read2: Path,
    targets: dict[int, set[int]],
) -> tuple[dict[int, Counter[int]], int]:
    support: dict[int, Counter[int]] = {k: Counter() for k in targets}
    total_pairs = 0
    for left, right in zip(fastq_records(read1), fastq_records(read2)):
        total_pairs += 1
        for k, wanted in targets.items():
            fragment = rolling_target_keys(left[0], k, wanted)
            fragment.update(rolling_target_keys(right[0], k, wanted))
            support[k].update(fragment)
    return support, total_pairs


def longest_run(ordered: list[int | None], support: Counter[int], threshold: int) -> int:
    best = 0
    current = 0
    for key in ordered:
        if key is not None and support.get(key, 0) >= threshold:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def summarize_seed(item: SeedKmers, support: Counter[int], membership: Counter[int]) -> dict[str, object]:
    keys = item.unique
    zero = sum(support.get(key, 0) == 0 for key in keys)
    singleton = sum(support.get(key, 0) == 1 for key in keys)
    solid = sum(support.get(key, 0) >= 2 for key in keys)
    shared = sum(membership.get(key, 0) > 1 for key in keys)
    n = max(1, len(keys))
    return {
        "seed_id": item.seed_id,
        "name": item.name,
        "length": item.length,
        "k": item.k,
        "distinct_target_kmers": len(keys),
        "zero_fragment_kmers": zero,
        "singleton_fragment_kmers": singleton,
        "solid_fragment_kmers": solid,
        "zero_fraction": zero / n,
        "singleton_fraction": singleton / n,
        "solid_fraction": solid / n,
        "shared_between_seed_kmers": shared,
        "shared_fraction": shared / n,
        "longest_raw_supported_run_kmers": longest_run(item.ordered, support, 1),
        "longest_solid_run_kmers": longest_run(item.ordered, support, 2),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--ks", default="31,41,55")
    ap.add_argument("--output-prefix", type=Path, required=True)
    args = ap.parse_args()

    ks = sorted({int(value) for value in args.ks.split(",") if value.strip()})
    seed_kmers, targets, membership = build_seed_kmers(args.seeds, ks)
    support, total_pairs = scan_fragment_support(args.read1, args.read2, targets)
    rows = [summarize_seed(item, support[item.k], membership[item.k]) for item in seed_kmers]

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    details = args.output_prefix.with_name(args.output_prefix.name + "_seeds.tsv")
    fields = list(rows[0]) if rows else []
    with details.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in ("zero_fraction", "singleton_fraction", "solid_fraction", "shared_fraction"):
                out[key] = f"{float(out[key]):.6f}"
            writer.writerow(out)

    summary_rows: list[dict[str, object]] = []
    for k in ks:
        wanted = targets[k]
        counts = support[k]
        zero = sum(counts.get(key, 0) == 0 for key in wanted)
        singleton = sum(counts.get(key, 0) == 1 for key in wanted)
        solid = sum(counts.get(key, 0) >= 2 for key in wanted)
        shared = sum(membership[k].get(key, 0) > 1 for key in wanted)
        per_seed = [row for row in rows if row["k"] == k]
        summary_rows.append(
            {
                "k": k,
                "seed_records": len(per_seed),
                "target_kmers": len(wanted),
                "zero_fragment_kmers": zero,
                "singleton_fragment_kmers": singleton,
                "solid_fragment_kmers": solid,
                "zero_fraction": zero / max(1, len(wanted)),
                "singleton_fraction": singleton / max(1, len(wanted)),
                "solid_fraction": solid / max(1, len(wanted)),
                "shared_between_seed_kmers": shared,
                "seeds_with_singleton_fraction_ge_0_10": sum(float(row["singleton_fraction"]) >= 0.10 for row in per_seed),
                "seeds_with_raw_supported_run_ge_50": sum(int(row["longest_raw_supported_run_kmers"]) >= 50 for row in per_seed),
                "seeds_with_solid_run_ge_50": sum(int(row["longest_solid_run_kmers"]) >= 50 for row in per_seed),
            }
        )

    summary_path = args.output_prefix.with_name(args.output_prefix.name + "_summary.tsv")
    sfields = list(summary_rows[0]) if summary_rows else []
    with summary_path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=sfields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            out = dict(row)
            for key in ("zero_fraction", "singleton_fraction", "solid_fraction"):
                out[key] = f"{float(out[key]):.6f}"
            writer.writerow(out)

    stats = {
        "read_pairs": total_pairs,
        "seed_records": len({item.seed_id for item in seed_kmers}),
        "ks": ks,
        "summary": summary_rows,
    }
    json_path = args.output_prefix.with_name(args.output_prefix.name + "_stats.json")
    json_path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(summary_path.read_text())


if __name__ == "__main__":
    main()
