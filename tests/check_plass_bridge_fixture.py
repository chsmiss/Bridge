#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1].upper()


def read_fasta(path: Path) -> List[str]:
    sequences: List[str] = []
    chunks: List[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(">"):
            if chunks:
                sequences.append("".join(chunks).upper())
            chunks = []
        else:
            chunks.append(raw_line.strip())
    if chunks:
        sequences.append("".join(chunks).upper())
    return sequences


def read_tsv(path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            (row["source"], row["target"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--plass-proteins", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    expected_records = read_fasta(args.expected)
    contigs = read_fasta(args.contigs)
    plass_proteins = read_fasta(args.plass_proteins)
    if len(expected_records) != 1:
        raise SystemExit("expected FASTA must contain exactly one sequence")
    if not contigs:
        raise SystemExit("no evidence-path contigs were produced")
    if not plass_proteins:
        raise SystemExit("Plass produced no protein sequences")

    expected = expected_records[0]
    canonical_expected = min(expected, reverse_complement(expected))
    if canonical_expected not in contigs:
        raise SystemExit(
            f"true reconstructed path not found; expected length={len(expected)} "
            f"observed lengths={[len(sequence) for sequence in contigs[:10]]}"
        )

    evidence = read_tsv(args.evidence)
    report = read_tsv(args.report)
    true_row = evidence.get(("A", "B"))
    decoy_row = evidence.get(("A", "C"))
    if true_row is None or decoy_row is None:
        raise SystemExit("missing expected A->B or A->C evidence rows")
    if true_row["breakpoint_class"] != "same_orf_supported":
        raise SystemExit(f"true edge was not protein-supported: {true_row}")
    if float(true_row["protein_score"]) <= float(decoy_row["protein_score"]):
        raise SystemExit("true edge protein score did not exceed decoy score")
    if report.get(("A", "B"), {}).get("selected") != "true":
        raise SystemExit("evidence path did not select A->B")
    if report.get(("A", "C"), {}).get("selected") != "false":
        raise SystemExit("evidence path selected the decoy A->C edge")

    print(
        "Plass bridge fixture passed: "
        f"proteins={len(plass_proteins)} true_score={true_row['protein_score']} "
        f"decoy_class={decoy_row['breakpoint_class']} contigs={len(contigs)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
