#!/usr/bin/env python3
"""Predict a sparse ORF set for BridgeBin's protein Biological Brain.

This is intentionally a *candidate-gated* gene caller.  When ``--pairs`` is supplied only
contigs that occur as candidate endpoints are processed.  In each contig we keep the
longest translated ORFs, because ESM-C is reserved for difficult bins rather than run on
every predicted protein in the metagenome.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import pyrodigal


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs", type=Path, required=True)
    p.add_argument("--output-faa", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--pairs", type=Path, help="optional candidate table; only pair endpoints are called")
    p.add_argument("--min-aa", type=int, default=50)
    p.add_argument("--max-proteins-per-contig", type=int, default=6)
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


def endpoints(path: Optional[Path]) -> Optional[Set[str]]:
    if path is None:
        return None
    selected: Set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for key in ("left", "right", "source", "target", "contig_a", "contig_b"):
                value = (row.get(key) or "").strip()
                if value:
                    selected.add(value)
    return selected


def clean_translation(value: object) -> str:
    seq = str(value).replace("*", "").upper()
    return "".join(residue if residue.isalpha() else "X" for residue in seq)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.min_aa < 10 or args.max_proteins_per_contig < 1:
        raise SystemExit("--min-aa must be >=10 and --max-proteins-per-contig positive")
    selected = endpoints(args.pairs)
    finder = pyrodigal.GeneFinder(meta=True)
    args.output_faa.parent.mkdir(parents=True, exist_ok=True)
    args.mapping.parent.mkdir(parents=True, exist_ok=True)
    contigs = proteins = 0
    with args.output_faa.open("w", encoding="utf-8") as faa, args.mapping.open(
        "w", encoding="utf-8", newline=""
    ) as mapping:
        writer = csv.writer(mapping, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "protein_id", "aa_length", "rank"])
        for contig, sequence in fasta(args.contigs):
            if selected is not None and contig not in selected:
                continue
            contigs += 1
            predicted = []
            for index, gene in enumerate(finder.find_genes(sequence.encode("ascii", errors="ignore")), start=1):
                protein = clean_translation(gene.translate(include_stop=False))
                if len(protein) >= args.min_aa:
                    predicted.append((len(protein), index, protein))
            predicted.sort(key=lambda item: (-item[0], item[1]))
            for rank, (aa_length, original_index, protein) in enumerate(
                predicted[: args.max_proteins_per_contig], start=1
            ):
                protein_id = f"{contig}::orf{original_index:04d}"
                faa.write(f">{protein_id}\n")
                for start in range(0, len(protein), 80):
                    faa.write(protein[start : start + 80] + "\n")
                writer.writerow([contig, protein_id, aa_length, rank])
                proteins += 1
    print(f"bridgebin-orfs: contigs={contigs} proteins={proteins} output={args.output_faa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
