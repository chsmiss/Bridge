#!/usr/bin/env python3
"""Benchmark-only oracle profiler for Stage24 one-fragment singleton rescue.

This script never participates in production assembly. It maps Stage24 guided
segments to known benchmark references and asks which reference-free segment
features enrich for reference-consistent one-fragment rescues.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional

BACTERIA = (
    "Bacillus_subtilis",
    "Enterococcus_faecalis",
    "Escherichia_coli",
    "Lactobacillus_fermentum",
    "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa",
    "Salmonella_enterica",
    "Staphylococcus_aureus",
)

HEADER_RE = re.compile(
    r"source=(\S+):(\d+)-(\d+)\s+"
    r"singleton_keys=(\d+)\s+singleton_fragments=(\d+)\s+solid_nodes=(\d+)"
)


@dataclass
class Segment:
    name: str
    sequence: str
    source: str
    source_start: int
    source_end: int
    singleton_keys: int
    singleton_fragments: int
    solid_nodes: int

    @property
    def length(self) -> int:
        return len(self.sequence)

    def singleton_density(self, k: int) -> float:
        return self.singleton_keys / max(1, self.length - k + 1)

    @property
    def singleton_fraction(self) -> float:
        return self.singleton_keys / max(1, self.singleton_keys + self.solid_nodes)


@dataclass
class Alignment:
    query: str
    qlen: int
    qstart: int
    qend: int
    target: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alen: int
    mapq: int
    species: str

    @property
    def query_coverage(self) -> float:
        return (self.qend - self.qstart) / max(1, self.qlen)

    @property
    def identity(self) -> float:
        return self.nmatch / max(1, self.alen)

    @property
    def score(self) -> tuple[float, float, int, int]:
        return (
            self.query_coverage * self.identity,
            self.query_coverage,
            self.mapq,
            self.nmatch,
        )


def read_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    name: Optional[str] = None
    header = ""
    seq_parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(seq_parts).upper()
                header = line[1:]
                name = header.split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        yield name, header, "".join(seq_parts).upper()


def parse_segments(path: Path) -> Dict[str, Segment]:
    out: Dict[str, Segment] = {}
    for name, header, seq in read_fasta(path):
        match = HEADER_RE.search(header)
        if not match:
            raise ValueError(f"unrecognized Stage24 segment header: {header}")
        source, start, end, keys, fragments, solid = match.groups()
        out[name] = Segment(
            name=name,
            sequence=seq,
            source=source,
            source_start=int(start),
            source_end=int(end),
            singleton_keys=int(keys),
            singleton_fragments=int(fragments),
            solid_nodes=int(solid),
        )
    return out


def parse_paf(path: Path) -> Dict[str, Alignment]:
    best: Dict[str, Alignment] = {}
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            target = fields[5]
            species = target.split("|", 1)[0]
            aln = Alignment(
                query=fields[0],
                qlen=int(fields[1]),
                qstart=int(fields[2]),
                qend=int(fields[3]),
                target=target,
                tlen=int(fields[6]),
                tstart=int(fields[7]),
                tend=int(fields[8]),
                nmatch=int(fields[9]),
                alen=int(fields[10]),
                mapq=int(fields[11]),
                species=species,
            )
            previous = best.get(aln.query)
            if previous is None or aln.score > previous.score:
                best[aln.query] = aln
    return best


def combine_refs(ref_dir: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with output.open("w") as out:
        for species in BACTERIA:
            path = ref_dir / f"{species}.fasta"
            if not path.exists():
                raise FileNotFoundError(path)
            for name, _header, seq in read_fasta(path):
                out.write(f">{species}|{name}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + "\n")
                written += 1
    if written != len(BACTERIA):
        raise RuntimeError(f"expected {len(BACTERIA)} bacterial references, wrote {written}")


def quantile(values: Iterable[float], p: float) -> float:
    xs = sorted(values)
    if not xs:
        return 0.0
    pos = (len(xs) - 1) * p
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return float(xs[lo])
    return float(xs[lo] * (hi - pos) + xs[hi] * (pos - lo))


def truth_flags(
    aln: Optional[Alignment],
    min_qcov: float,
    min_identity: float,
    strict_qcov: float,
    strict_identity: float,
) -> tuple[bool, bool]:
    if aln is None:
        return False, False
    consistent = aln.query_coverage >= min_qcov and aln.identity >= min_identity
    strict = aln.query_coverage >= strict_qcov and aln.identity >= strict_identity
    return consistent, strict


def rule_scan(rows: list[dict], output: Path) -> list[dict]:
    one = [r for r in rows if r["singleton_fragments"] == 1]
    total_true = sum(r["reference_consistent"] for r in one)
    densities = (0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 1.01)
    min_solids = (2, 20, 40, 60, 80, 100, 150)
    max_keys = (5, 10, 20, 30, 40, 60, 10**9)
    rules: list[dict] = []
    for density in densities:
        for min_solid in min_solids:
            for max_key in max_keys:
                selected = [
                    r
                    for r in one
                    if r["singleton_density"] <= density
                    and r["solid_nodes"] >= min_solid
                    and r["singleton_keys"] <= max_key
                ]
                if not selected:
                    continue
                true = sum(r["reference_consistent"] for r in selected)
                strict = sum(r["strict_reference_consistent"] for r in selected)
                rules.append(
                    {
                        "max_singleton_density": density,
                        "min_solid_nodes": min_solid,
                        "max_singleton_keys": max_key if max_key < 10**9 else -1,
                        "selected": len(selected),
                        "selected_bases": sum(r["length"] for r in selected),
                        "reference_consistent": true,
                        "strict_reference_consistent": strict,
                        "precision": true / len(selected),
                        "recall": true / total_true if total_true else 0.0,
                    }
                )
    rules.sort(key=lambda r: (r["precision"], r["reference_consistent"], r["selected"]), reverse=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        fieldnames = list(rules[0].keys()) if rules else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if fieldnames:
            writer.writeheader()
            writer.writerows(rules)
    return rules


def analyze(args: argparse.Namespace) -> None:
    segments = parse_segments(args.segments)
    alignments = parse_paf(args.paf)
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for name, segment in segments.items():
        aln = alignments.get(name)
        consistent, strict = truth_flags(
            aln,
            args.min_query_coverage,
            args.min_identity,
            args.strict_query_coverage,
            args.strict_identity,
        )
        row = {
            "name": name,
            "length": segment.length,
            "source": segment.source,
            "singleton_keys": segment.singleton_keys,
            "singleton_fragments": segment.singleton_fragments,
            "solid_nodes": segment.solid_nodes,
            "singleton_density": segment.singleton_density(args.k),
            "singleton_fraction": segment.singleton_fraction,
            "mapped": aln is not None,
            "species": aln.species if aln else "",
            "query_coverage": aln.query_coverage if aln else 0.0,
            "identity": aln.identity if aln else 0.0,
            "mapq": aln.mapq if aln else 0,
            "reference_consistent": consistent,
            "strict_reference_consistent": strict,
        }
        rows.append(row)

    fields = list(rows[0].keys()) if rows else []
    with (outdir / "segment_oracle.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        if fields:
            writer.writeheader()
            writer.writerows(rows)

    summary: dict[str, object] = {
        "policy": {
            "benchmark_only": True,
            "production_reference_free": True,
            "k": args.k,
            "reference_consistent": {
                "min_query_coverage": args.min_query_coverage,
                "min_identity": args.min_identity,
            },
            "strict_reference_consistent": {
                "min_query_coverage": args.strict_query_coverage,
                "min_identity": args.strict_identity,
            },
        },
        "groups": {},
        "one_fragment_species": {},
    }

    for group_name, predicate in (
        ("one_fragment", lambda r: r["singleton_fragments"] == 1),
        ("multi_fragment", lambda r: r["singleton_fragments"] >= 2),
        ("all", lambda r: True),
    ):
        group = [r for r in rows if predicate(r)]
        true = [r for r in group if r["reference_consistent"]]
        strict = [r for r in group if r["strict_reference_consistent"]]
        mapped = [r for r in group if r["mapped"]]
        summary["groups"][group_name] = {
            "segments": len(group),
            "bases": sum(r["length"] for r in group),
            "mapped": len(mapped),
            "reference_consistent": len(true),
            "strict_reference_consistent": len(strict),
            "reference_consistent_rate": len(true) / len(group) if group else 0.0,
            "strict_rate": len(strict) / len(group) if group else 0.0,
            "length_median": quantile((r["length"] for r in group), 0.5),
            "length_p90": quantile((r["length"] for r in group), 0.9),
            "singleton_density_median": quantile((r["singleton_density"] for r in group), 0.5),
            "true_length_median": quantile((r["length"] for r in true), 0.5),
            "true_singleton_density_median": quantile((r["singleton_density"] for r in true), 0.5),
        }

    species = Counter(
        r["species"]
        for r in rows
        if r["singleton_fragments"] == 1 and r["reference_consistent"] and r["species"]
    )
    summary["one_fragment_species"] = dict(species.most_common())

    rules = rule_scan(rows, outdir / "one_fragment_rule_scan.tsv")
    shortlist = []
    for min_recall in (0.20, 0.40, 0.60, 0.80):
        candidates = [r for r in rules if r["recall"] >= min_recall]
        if candidates:
            best = max(candidates, key=lambda r: (r["precision"], r["reference_consistent"], -r["selected"]))
            item = dict(best)
            item["minimum_recall_bucket"] = min_recall
            shortlist.append(item)
    summary["rule_shortlist"] = shortlist

    with (outdir / "stage25_oracle_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    with (outdir / "stage25_oracle_summary.tsv").open("w") as handle:
        handle.write("group\tsegments\tbases\tmapped\treference_consistent\tstrict_reference_consistent\treference_consistent_rate\tstrict_rate\tlength_median\tlength_p90\tsingleton_density_median\ttrue_length_median\ttrue_singleton_density_median\n")
        for group_name in ("one_fragment", "multi_fragment", "all"):
            g = summary["groups"][group_name]
            handle.write(
                f"{group_name}\t{g['segments']}\t{g['bases']}\t{g['mapped']}\t{g['reference_consistent']}\t{g['strict_reference_consistent']}\t"
                f"{g['reference_consistent_rate']:.6f}\t{g['strict_rate']:.6f}\t{g['length_median']:.1f}\t{g['length_p90']:.1f}\t"
                f"{g['singleton_density_median']:.6f}\t{g['true_length_median']:.1f}\t{g['true_singleton_density_median']:.6f}\n"
            )

    print((outdir / "stage25_oracle_summary.tsv").read_text(), end="")
    print("one_fragment_species", json.dumps(summary["one_fragment_species"], sort_keys=True))
    print("rule_shortlist", json.dumps(shortlist, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    combine = sub.add_parser("combine-refs")
    combine.add_argument("--ref-dir", type=Path, required=True)
    combine.add_argument("--output", type=Path, required=True)

    profile = sub.add_parser("analyze")
    profile.add_argument("--segments", type=Path, required=True)
    profile.add_argument("--paf", type=Path, required=True)
    profile.add_argument("--output-dir", type=Path, required=True)
    profile.add_argument("--k", type=int, default=31)
    profile.add_argument("--min-query-coverage", type=float, default=0.80)
    profile.add_argument("--min-identity", type=float, default=0.97)
    profile.add_argument("--strict-query-coverage", type=float, default=0.90)
    profile.add_argument("--strict-identity", type=float, default=0.99)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "combine-refs":
        combine_refs(args.ref_dir, args.output)
    elif args.command == "analyze":
        analyze(args)
    else:
        raise AssertionError(args.command)


if __name__ == "__main__":
    main()
