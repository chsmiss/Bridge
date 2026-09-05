#!/usr/bin/env python3
"""Summarize reference-protein recovery from MMseqs2 tabular alignments.

The expected alignment format is produced by::

    mmseqs easy-search reference.faa predicted.faa hits.m8 tmp \
      --format-output query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits

Two recall views are reported:

* single-hit recall: one predicted protein spans the requested fraction of a
  reference protein;
* union recall: multiple fragments may jointly cover the reference protein.

Reciprocal completeness additionally requires target coverage, making it a
more conservative proxy for complete protein recovery.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Hit:
    query: str
    target: str
    identity: float
    query_start: int
    query_end: int
    query_length: int
    target_start: int
    target_end: int
    target_length: int

    @property
    def query_coverage(self) -> float:
        return (self.query_end - self.query_start) / self.query_length

    @property
    def target_coverage(self) -> float:
        return (self.target_end - self.target_start) / self.target_length


@dataclass
class SampleStats:
    label: str
    predicted_proteins: int
    predicted_amino_acids: int
    predicted_ge_100aa: int
    alignment_records: int
    identity_passing_hits: int
    reference_proteins: int
    reference_amino_acids: int
    recalled_single: int
    recalled_union: int
    recalled_reciprocal: int
    reference_covered_amino_acids: int
    matched_predicted: int
    single_recall_fraction: float
    union_recall_fraction: float
    reciprocal_recall_fraction: float
    reference_aa_coverage_fraction: float
    matched_predicted_fraction: float


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    current = 0
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = current
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty FASTA header in {path} at line {line_number}")
                name = header.split()[0]
                if name in lengths:
                    raise ValueError(f"duplicate FASTA identifier {name!r} in {path}")
                current = 0
            else:
                if name is None:
                    raise ValueError(
                        f"sequence before first FASTA header in {path} at line {line_number}"
                    )
                current += len(line)
    if name is not None:
        lengths[name] = current
    if not lengths:
        raise ValueError(f"no FASTA records found in {path}")
    empty = [record for record, length in lengths.items() if length == 0]
    if empty:
        raise ValueError(f"empty FASTA records in {path}: {empty[:5]}")
    return lengths


def normalized_identity(value: str) -> float:
    identity = float(value)
    if identity > 1.0:
        identity /= 100.0
    if not 0.0 <= identity <= 1.0:
        raise ValueError(f"identity outside 0..1 after normalization: {value}")
    return identity


def zero_based_interval(start_text: str, end_text: str, length: int) -> tuple[int, int]:
    start = int(start_text)
    end = int(end_text)
    low = min(start, end)
    high = max(start, end)
    # MMseqs2 reports 1-based inclusive coordinates. The fallback accepts a
    # zero start as a zero-based coordinate for synthetic tests.
    if low > 0:
        low -= 1
    high = min(high, length)
    low = max(0, min(low, length))
    if high <= low:
        raise ValueError(f"invalid alignment interval {start_text}:{end_text}/{length}")
    return low, high


def read_hits(
    path: Path,
    references: dict[str, int],
    predictions: dict[str, int],
) -> tuple[list[Hit], int]:
    hits: list[Hit] = []
    records = 0
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            records += 1
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"{path}:{line_number} has {len(fields)} fields; expected at least 12"
                )
            query, target = fields[0], fields[1]
            if query not in references or target not in predictions:
                continue
            query_length = references[query]
            target_length = predictions[target]
            reported_query_length = int(fields[6])
            reported_target_length = int(fields[9])
            if reported_query_length != query_length:
                raise ValueError(
                    f"query length mismatch for {query}: FASTA={query_length}, "
                    f"alignment={reported_query_length}"
                )
            if reported_target_length != target_length:
                raise ValueError(
                    f"target length mismatch for {target}: FASTA={target_length}, "
                    f"alignment={reported_target_length}"
                )
            query_start, query_end = zero_based_interval(
                fields[4], fields[5], query_length
            )
            target_start, target_end = zero_based_interval(
                fields[7], fields[8], target_length
            )
            hits.append(
                Hit(
                    query=query,
                    target=target,
                    identity=normalized_identity(fields[2]),
                    query_start=query_start,
                    query_end=query_end,
                    query_length=query_length,
                    target_start=target_start,
                    target_end=target_end,
                    target_length=target_length,
                )
            )
    return hits, records


def merged_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def safe_fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def summarize_sample(
    label: str,
    predicted_fasta: Path,
    alignment_path: Path,
    references: dict[str, int],
    min_identity: float,
    min_query_coverage: float,
    min_target_coverage: float,
) -> SampleStats:
    predictions = fasta_lengths(predicted_fasta)
    hits, records = read_hits(alignment_path, references, predictions)
    passing = [hit for hit in hits if hit.identity >= min_identity]

    by_query: dict[str, list[Hit]] = defaultdict(list)
    matched_targets: set[str] = set()
    for hit in passing:
        by_query[hit.query].append(hit)
        if hit.target_coverage >= min_target_coverage:
            matched_targets.add(hit.target)

    recalled_single = 0
    recalled_union = 0
    recalled_reciprocal = 0
    covered_amino_acids = 0
    for reference, reference_length in references.items():
        query_hits = by_query.get(reference, [])
        if any(hit.query_coverage >= min_query_coverage for hit in query_hits):
            recalled_single += 1
        if any(
            hit.query_coverage >= min_query_coverage
            and hit.target_coverage >= min_target_coverage
            for hit in query_hits
        ):
            recalled_reciprocal += 1
        union_length = merged_length(
            [(hit.query_start, hit.query_end) for hit in query_hits]
        )
        covered_amino_acids += union_length
        if union_length / reference_length >= min_query_coverage:
            recalled_union += 1

    reference_count = len(references)
    reference_aa = sum(references.values())
    predicted_count = len(predictions)
    return SampleStats(
        label=label,
        predicted_proteins=predicted_count,
        predicted_amino_acids=sum(predictions.values()),
        predicted_ge_100aa=sum(length >= 100 for length in predictions.values()),
        alignment_records=records,
        identity_passing_hits=len(passing),
        reference_proteins=reference_count,
        reference_amino_acids=reference_aa,
        recalled_single=recalled_single,
        recalled_union=recalled_union,
        recalled_reciprocal=recalled_reciprocal,
        reference_covered_amino_acids=covered_amino_acids,
        matched_predicted=len(matched_targets),
        single_recall_fraction=safe_fraction(recalled_single, reference_count),
        union_recall_fraction=safe_fraction(recalled_union, reference_count),
        reciprocal_recall_fraction=safe_fraction(
            recalled_reciprocal, reference_count
        ),
        reference_aa_coverage_fraction=safe_fraction(
            covered_amino_acids, reference_aa
        ),
        matched_predicted_fraction=safe_fraction(
            len(matched_targets), predicted_count
        ),
    )


def print_tsv(stats: list[SampleStats]) -> None:
    fields = list(SampleStats.__dataclass_fields__)
    print("\t".join(fields))
    for sample in stats:
        values: list[str] = []
        for field in fields:
            value = getattr(sample, field)
            values.append(f"{value:.6f}" if isinstance(value, float) else str(value))
        print("\t".join(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument(
        "--sample",
        nargs=3,
        action="append",
        metavar=("LABEL", "PREDICTED_FASTA", "ALIGNMENTS_M8"),
        required=True,
    )
    parser.add_argument("--min-identity", type=float, default=0.80)
    parser.add_argument("--min-query-coverage", type=float, default=0.80)
    parser.add_argument("--min-target-coverage", type=float, default=0.80)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    for name, value in (
        ("--min-identity", args.min_identity),
        ("--min-query-coverage", args.min_query_coverage),
        ("--min-target-coverage", args.min_target_coverage),
    ):
        if not 0.0 <= value <= 1.0:
            parser.error(f"{name} must be in 0..1")

    references = fasta_lengths(args.reference)
    labels: set[str] = set()
    output: list[SampleStats] = []
    for label, predicted_text, alignment_text in args.sample:
        if label in labels:
            parser.error(f"duplicate sample label: {label}")
        labels.add(label)
        output.append(
            summarize_sample(
                label=label,
                predicted_fasta=Path(predicted_text),
                alignment_path=Path(alignment_text),
                references=references,
                min_identity=args.min_identity,
                min_query_coverage=args.min_query_coverage,
                min_target_coverage=args.min_target_coverage,
            )
        )

    print_tsv(output)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(
            json.dumps([asdict(sample) for sample in output], indent=2, sort_keys=True)
            + "\n"
        )


if __name__ == "__main__":
    main()
