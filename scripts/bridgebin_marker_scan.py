#!/usr/bin/env python3
"""Single-copy-marker contamination focus for BridgeBin.

This module deliberately does *not* assign genomes.  It answers three narrower questions:

1. Which contigs contain trusted single-copy marker hits?
2. Which contig pairs in the same current bin are biologically incompatible because they
   carry the same single-copy marker?
3. Which current bins deserve expensive DNA/gene/protein refinement, and what is a
   conservative lower bound on the number of genomes represented by marker multiplicity?

The default marker HMM is the 107-marker collection distributed by SemiBin2 (MIT).  It is
not vendored into Bridge because it is ~16 MB; this script can download the current file
or consume an explicit ``--marker-hmm`` path.  ORFs are predicted with Pyrodigal and all
ORFs are scanned with HMMER ``hmmsearch --cut_tc``.  Marker domains covering <=40% of the
HMM are rejected, matching the conservative filtering used by SemiBin2.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import subprocess
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import pyrodigal

SEMIBIN_MARKER_URL = (
    "https://raw.githubusercontent.com/BigDataBiology/SemiBin/main/SemiBin/marker.hmm"
)
NORMALIZE_MARKER = {
    "TIGR00388": "TIGR00389",
    "TIGR00471": "TIGR00472",
    "TIGR00408": "TIGR00409",
    "TIGR02386": "TIGR02387",
}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--assignments", type=Path, help="optional current assignments.tsv")
    p.add_argument("--marker-hmm", type=Path, help="local marker HMM; otherwise download SemiBin2 marker.hmm")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-contig", type=int, default=1500)
    p.add_argument("--min-marker-coverage", type=float, default=0.40)
    p.add_argument("--focus-min-duplicate-markers", type=int, default=1)
    p.add_argument("--focus-min-duplicate-fraction", type=float, default=0.025)
    p.add_argument("--keep-proteins", action="store_true")
    return p.parse_args(argv)


def fasta(path: Path) -> Iterator[Tuple[str, str]]:
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
                    raise ValueError(f"{path}: sequence before FASTA header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def read_assignments(path: Optional[Path]) -> Dict[str, Optional[str]]:
    if path is None:
        return {}
    result: Dict[str, Optional[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "contig" not in reader.fieldnames or "bin" not in reader.fieldnames:
            raise ValueError(f"{path}: need contig and bin columns")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            raw_bin = (row.get("bin") or "").strip()
            if contig:
                result[contig] = None if raw_bin in {"", ".", "NA", "unbinned"} else raw_bin
    return result


def write_all_orfs(contigs: Path, proteins: Path, mapping: Path, min_contig: int) -> Tuple[int, int]:
    finder = pyrodigal.GeneFinder(meta=True)
    proteins.parent.mkdir(parents=True, exist_ok=True)
    contig_count = protein_count = 0
    with proteins.open("w", encoding="utf-8") as faa, mapping.open(
        "w", encoding="utf-8", newline=""
    ) as map_handle:
        writer = csv.writer(map_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "contig", "aa_length"])
        for contig, sequence in fasta(contigs):
            if len(sequence) < min_contig:
                continue
            contig_count += 1
            genes = finder.find_genes(sequence.encode("ascii", errors="ignore"))
            for index, gene in enumerate(genes, start=1):
                protein = str(gene.translate(include_stop=False)).rstrip("*").upper()
                if not protein:
                    continue
                protein_id = f"bb_orf_{protein_count:09d}"
                faa.write(f">{protein_id}\n")
                for start in range(0, len(protein), 80):
                    faa.write(protein[start : start + 80] + "\n")
                writer.writerow([protein_id, contig, len(protein)])
                protein_count += 1
    return contig_count, protein_count


def ensure_marker_hmm(path: Optional[Path], output_dir: Path) -> Path:
    if path is not None:
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    target = output_dir / "semibin107.marker.hmm"
    if not target.exists() or target.stat().st_size < 1000000:
        print(f"bridgebin-markers: downloading marker HMM from {SEMIBIN_MARKER_URL}")
        urllib.request.urlretrieve(SEMIBIN_MARKER_URL, target)
    return target


def read_orf_mapping(path: Path) -> Dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["protein_id"]: row["contig"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("protein_id") and row.get("contig")
        }


def parse_domtbl(
    path: Path, orf_to_contig: Dict[str, str], min_marker_coverage: float
) -> List[Dict[str, object]]:
    best: Dict[Tuple[str, str], Dict[str, object]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 22:
                continue
            protein_id = fields[0]
            contig = orf_to_contig.get(protein_id)
            if contig is None:
                continue
            marker = NORMALIZE_MARKER.get(fields[3], fields[3])
            try:
                qlen = int(fields[5])
                i_evalue = float(fields[12])
                score = float(fields[13])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
            except ValueError:
                continue
            coverage = max(0.0, (hmm_to - hmm_from + 1) / max(1, qlen))
            if coverage <= min_marker_coverage:
                continue
            row = {
                "contig": contig,
                "protein_id": protein_id,
                "marker": marker,
                "marker_coverage": coverage,
                "i_evalue": i_evalue,
                "domain_score": score,
            }
            key = (contig, marker)
            old = best.get(key)
            if old is None or score > float(old["domain_score"]):
                best[key] = row
    return sorted(best.values(), key=lambda row: (str(row["contig"]), str(row["marker"])))


def write_hits(path: Path, hits: List[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["contig", "protein_id", "marker", "marker_coverage", "i_evalue", "domain_score"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for hit in hits:
            row = dict(hit)
            row["marker_coverage"] = f"{float(row['marker_coverage']):.6f}"
            row["i_evalue"] = f"{float(row['i_evalue']):.6g}"
            row["domain_score"] = f"{float(row['domain_score']):.3f}"
            writer.writerow(row)


def summarize(
    hits: List[Dict[str, object]],
    assignments: Dict[str, Optional[str]],
    report_path: Path,
    cannot_path: Path,
    focus_path: Path,
    focus_min_duplicate_markers: int,
    focus_min_duplicate_fraction: float,
) -> Dict[str, object]:
    marker_contigs: Dict[Tuple[str, str], set[str]] = defaultdict(set)
    observed_per_bin: Dict[str, Counter[str]] = defaultdict(Counter)
    hit_contigs_per_bin: Dict[str, set[str]] = defaultdict(set)
    unassigned_hits = 0
    for hit in hits:
        contig = str(hit["contig"])
        bin_name = assignments.get(contig) if assignments else "all"
        if bin_name is None:
            unassigned_hits += 1
            continue
        marker = str(hit["marker"])
        marker_contigs[(bin_name, marker)].add(contig)
        observed_per_bin[bin_name][marker] += 1
        hit_contigs_per_bin[bin_name].add(contig)

    cannot_rows = []
    for (bin_name, marker), contigs in sorted(marker_contigs.items()):
        if len(contigs) < 2:
            continue
        for left, right in itertools.combinations(sorted(contigs), 2):
            cannot_rows.append((left, right, bin_name, marker, "duplicate_single_copy_marker"))
    with cannot_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["left", "right", "bin", "marker", "reason"])
        writer.writerows(cannot_rows)

    bins = []
    focused = []
    for bin_name in sorted(observed_per_bin):
        counts = observed_per_bin[bin_name]
        observed = len(counts)
        duplicated = sum(value >= 2 for value in counts.values())
        duplicate_fraction = duplicated / observed if observed else 0.0
        max_copy = max(counts.values(), default=0)
        k_lower_bound = max(1, max_copy)
        marker_excess = sum(max(0, value - 1) for value in counts.values())
        focus = duplicated >= focus_min_duplicate_markers or duplicate_fraction >= focus_min_duplicate_fraction
        if focus:
            focused.append(bin_name)
        bins.append(
            {
                "bin": bin_name,
                "marker_observed": observed,
                "marker_hits": sum(counts.values()),
                "marker_contigs": len(hit_contigs_per_bin[bin_name]),
                "duplicate_markers": duplicated,
                "duplicate_fraction": duplicate_fraction,
                "marker_excess_copies": marker_excess,
                "max_marker_copy": max_copy,
                "genome_count_lower_bound": k_lower_bound,
                "focus": focus,
            }
        )
    report = {
        "marker_hits": len(hits),
        "bins_with_markers": len(bins),
        "focused_bins": focused,
        "cannot_links": len(cannot_rows),
        "unassigned_marker_hits": unassigned_hits,
        "bins": bins,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with focus_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["bin"])
        for bin_name in focused:
            writer.writerow([bin_name])
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.threads < 1 or args.min_contig < 1:
        raise SystemExit("--threads and --min-contig must be positive")
    if not 0.0 < args.min_marker_coverage <= 1.0:
        raise SystemExit("--min-marker-coverage must be in (0,1]")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    proteins = args.output_dir / "all_orfs.faa"
    mapping = args.output_dir / "orf_to_contig.tsv"
    domtbl = args.output_dir / "markers.domtblout"
    hmm = ensure_marker_hmm(args.marker_hmm, args.output_dir)
    contigs, proteins_count = write_all_orfs(args.contigs, proteins, mapping, args.min_contig)
    if proteins_count == 0:
        raise SystemExit("no ORFs predicted")
    log_path = args.output_dir / "hmmsearch.log"
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(
            [
                "hmmsearch",
                "--domtblout",
                str(domtbl),
                "--cut_tc",
                "--cpu",
                str(args.threads),
                str(hmm),
                str(proteins),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    hits = parse_domtbl(domtbl, read_orf_mapping(mapping), args.min_marker_coverage)
    write_hits(args.output_dir / "marker_hits.tsv", hits)
    assignments = read_assignments(args.assignments)
    report = summarize(
        hits,
        assignments,
        args.output_dir / "marker_report.json",
        args.output_dir / "marker_cannot_links.tsv",
        args.output_dir / "marker_focus_bins.tsv",
        args.focus_min_duplicate_markers,
        args.focus_min_duplicate_fraction,
    )
    if not args.keep_proteins:
        proteins.unlink(missing_ok=True)
    print(
        f"bridgebin-markers: contigs={contigs} proteins={proteins_count} hits={len(hits)} "
        f"focused_bins={len(report['focused_bins'])} cannot_links={report['cannot_links']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
