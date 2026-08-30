#!/usr/bin/env python3
from __future__ import annotations

import csv
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
PROTEIN_EVIDENCE = ROOT / "biological_brain" / "protein_bridge_evidence.py"
DNA_LM = ROOT / "biological_brain" / "kmer_junction_lm.py"

CODONS: Dict[str, str] = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def back_translate(protein: str) -> str:
    return "".join(CODONS[amino_acid] for amino_acid in protein)


def random_protein(rng: random.Random, length: int) -> str:
    return "M" + "".join(rng.choice(AA_ALPHABET) for _ in range(length - 1))


def write_fixture(directory: Path, duplicate_reference: bool = False) -> Tuple[Path, Path, str, str]:
    rng = random.Random(43)
    protein = random_protein(rng, 220)
    dna = back_translate(protein)
    split = 390
    overlap = 31
    source = dna[:split]
    target = dna[split - overlap :]

    decoy_protein = random_protein(rng, (len(dna) - split + 2) // 3 + 12)
    decoy_suffix = back_translate(decoy_protein)[: len(dna) - split]
    decoy = dna[split - overlap : split] + decoy_suffix

    gfa = directory / "graph.gfa"
    gfa.write_text(
        "\n".join(
            [
                "H\tVN:Z:1.0",
                f"S\tA\t{source}\tKC:f:30.0",
                f"S\tB\t{target}\tKC:f:30.0",
                f"S\tC\t{decoy}\tKC:f:30.0",
                f"L\tA\t+\tB\t+\t{overlap}M\tDR:i:1\tGR:i:0\tPE:i:1",
                f"L\tA\t+\tC\t+\t{overlap}M\tDR:i:1\tGR:i:0\tPE:i:1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    proteins = directory / "proteins.faa"
    text = f">true_protein\n{protein}\n"
    if duplicate_reference:
        text += f">true_protein_copy\n{protein}\n"
    proteins.write_text(text, encoding="utf-8")
    return gfa, proteins, dna, decoy


def run(command: List[str], cwd: Path = ROOT) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def read_tsv(path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            (row["source"], row["target"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_fasta(path: Path) -> List[str]:
    sequences: List[str] = []
    chunks: List[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(">"):
            if chunks:
                sequences.append("".join(chunks))
            chunks = []
        else:
            chunks.append(raw_line.strip())
    if chunks:
        sequences.append("".join(chunks))
    return sequences


class BiologicalBrainTests(unittest.TestCase):
    def test_protein_assembly_selects_true_existing_branch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bridge-protein-") as temp:
            directory = Path(temp)
            gfa, proteins, expected_dna, _ = write_fixture(directory)
            evidence = directory / "protein_evidence.tsv"
            run(
                [
                    sys.executable,
                    str(PROTEIN_EVIDENCE),
                    "--gfa",
                    str(gfa),
                    "--proteins",
                    str(proteins),
                    "--output",
                    str(evidence),
                    "--junction-nt",
                    "450",
                ]
            )
            rows = read_tsv(evidence)
            self.assertEqual(set(rows), {("A", "B"), ("A", "C")})
            true_row = rows[("A", "B")]
            decoy_row = rows[("A", "C")]
            self.assertEqual(true_row["breakpoint_class"], "same_orf_supported")
            self.assertGreater(float(true_row["protein_score"]), 0.70)
            self.assertGreaterEqual(int(true_row["left_kmers"]), 2)
            self.assertGreaterEqual(int(true_row["right_kmers"]), 2)
            self.assertNotEqual(decoy_row["breakpoint_class"], "same_orf_supported")
            self.assertEqual(float(decoy_row["protein_score"]), 0.0)

            output = directory / "contigs.fasta"
            report = directory / "path_report.tsv"
            run(
                [
                    "cargo",
                    "run",
                    "--quiet",
                    "--bin",
                    "bridgeasm-evidence-path",
                    "--",
                    "--gfa",
                    str(gfa),
                    "--edge-evidence",
                    str(evidence),
                    "--output",
                    str(output),
                    "--report",
                    str(report),
                    "--mode",
                    "balanced",
                    "--min-length",
                    "100",
                ]
            )
            sequences = read_fasta(output)
            self.assertTrue(sequences)
            canonical_expected = min(expected_dna, reverse_complement(expected_dna))
            self.assertIn(canonical_expected, sequences)
            report_rows = read_tsv(report)
            self.assertEqual(report_rows[("A", "B")]["selected"], "true")
            self.assertEqual(report_rows[("A", "C")]["selected"], "false")

    def test_duplicate_protein_is_marked_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bridge-ambiguous-") as temp:
            directory = Path(temp)
            gfa, proteins, _, _ = write_fixture(directory, duplicate_reference=True)
            evidence = directory / "protein_evidence.tsv"
            run(
                [
                    sys.executable,
                    str(PROTEIN_EVIDENCE),
                    "--gfa",
                    str(gfa),
                    "--proteins",
                    str(proteins),
                    "--output",
                    str(evidence),
                ]
            )
            true_row = read_tsv(evidence)[("A", "B")]
            self.assertEqual(true_row["breakpoint_class"], "ambiguous_homology")
            self.assertEqual(float(true_row["protein_score"]), 0.0)
            self.assertGreater(float(true_row["ambiguity"]), 0.45)

    def test_dna_lm_adapter_emits_fusion_schema(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bridge-dnalm-") as temp:
            directory = Path(temp)
            gfa, _, _, _ = write_fixture(directory)
            output = directory / "dna_lm.tsv"
            run(
                [
                    sys.executable,
                    str(DNA_LM),
                    "--gfa",
                    str(gfa),
                    "--output",
                    str(output),
                ]
            )
            rows = read_tsv(output)
            self.assertEqual(set(rows), {("A", "B"), ("A", "C")})
            for row in rows.values():
                self.assertTrue(math_is_finite(row["dna_lm_delta"]))
                self.assertTrue(row["model"].startswith("markov_order_"))


def math_is_finite(value: str) -> bool:
    try:
        number = float(value)
    except ValueError:
        return False
    return number == number and number not in (float("inf"), float("-inf"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
