#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def run_script(name: str, *arguments: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(REPO / "scripts" / name), *(str(value) for value in arguments)],
        check=True,
        text=True,
        capture_output=True,
    )


class JunctionToolTests(unittest.TestCase):
    def test_plass_same_protein_support(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            backbone = root / "backbone.fasta"
            proteins = root / "backbone.faa"
            alignments = root / "backbone-to-plass.m8"
            edges = root / "edges.tsv"
            scores = root / "scores.tsv"
            details = root / "details.tsv"

            backbone.write_text(
                ">c1\n" + "A" * 300 + "\n>c2\n" + "C" * 300 + "\n"
            )
            proteins.write_text(
                ">c1_1 # 241 # 300 # 1 # ID=1_1\n"
                + "M" * 20
                + "\n>c2_1 # 1 # 60 # 1 # ID=2_1\n"
                + "M" * 20
                + "\n"
            )
            alignments.write_text(
                "c1_1\tprotein_x\t0.95\t20\t1\t20\t20\t1\t20\t40\t1e-20\t100\n"
                "c2_1\tprotein_x\t0.95\t20\t1\t20\t20\t21\t40\t40\t1e-20\t100\n"
            )
            edges.write_text(
                "source\ttarget\teligible\tprojected_gap\tguide_bases\n"
                "c1+\tc2+\ttrue\t0\t0\n"
            )

            result = run_script(
                "plass_junction_scores.py",
                "--backbone",
                backbone,
                "--proteins",
                proteins,
                "--alignments",
                alignments,
                "--edge-report",
                edges,
                "--output",
                scores,
                "--details",
                details,
            )
            rows = scores.read_text().strip().splitlines()
            self.assertEqual(len(rows), 2)
            fields = rows[1].split("\t")
            self.assertEqual(fields[0:2], ["c1+", "c2+"])
            self.assertGreater(float(fields[2]), 0.8)
            self.assertEqual(fields[4], "plass_same_protein")
            self.assertIn("supported_edges=1", result.stderr)

    def test_exact_penguin_bridge_window_export(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gfa = root / "backbone.gfa"
            guide = root / "guide.fasta"
            paf = root / "anchors.paf"
            edges = root / "edges.tsv"
            output = root / "junctions.fasta"
            manifest = root / "junctions.tsv"

            gfa.write_text(
                "H\tVN:Z:1.0\n"
                + "S\tc1\t"
                + "A" * 100
                + "\tLN:i:100\tKC:f:1.0\n"
                + "S\tc2\t"
                + "C" * 100
                + "\tLN:i:100\tKC:f:1.0\n"
            )
            guide.write_text(">g1\n" + "A" * 100 + "GGG" + "C" * 100 + "\n")
            paf.write_text(
                "c1\t100\t0\t100\t+\tg1\t203\t0\t100\t100\t100\t60\n"
                "c2\t100\t0\t100\t+\tg1\t203\t103\t203\t100\t100\t60\n"
            )
            edges.write_text(
                "source\ttarget\tguide\teligible\tselected\tprojected_gap\toverlap\tguide_bases\n"
                "c1+\tc2+\tg1\ttrue\ttrue\t3\t0\t3\n"
            )

            result = run_script(
                "export_junction_windows.py",
                "--gfa",
                gfa,
                "--guide",
                guide,
                "--paf",
                paf,
                "--edge-report",
                edges,
                "--output-fasta",
                output,
                "--manifest",
                manifest,
                "--left-context",
                10,
                "--right-context",
                10,
            )
            sequence = "".join(
                line.strip() for line in output.read_text().splitlines() if not line.startswith(">")
            )
            self.assertEqual(sequence, "A" * 10 + "GGG" + "C" * 10)
            manifest_row = manifest.read_text().strip().splitlines()[1].split("\t")
            header = manifest.read_text().strip().splitlines()[0].split("\t")
            record = dict(zip(header, manifest_row))
            self.assertEqual(record["bridge_bases"], "3")
            self.assertEqual(record["window_bases"], "23")
            self.assertIn("exported=1", result.stdout)

    def test_protein_recall_single_and_union_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.faa"
            predicted = root / "predicted.faa"
            alignments = root / "matches.m8"
            output_json = root / "recall.json"

            reference.write_text(">ref1\n" + "M" * 100 + "\n")
            predicted.write_text(">pred1\n" + "M" * 100 + "\n")
            alignments.write_text(
                "ref1\tpred1\t0.90\t90\t1\t90\t100\t1\t90\t100\t1e-20\t100\n"
            )
            run_script(
                "protein_recall.py",
                "--reference",
                reference,
                "--sample",
                "sample",
                predicted,
                alignments,
                "--min-identity",
                0.8,
                "--min-query-coverage",
                0.8,
                "--min-target-coverage",
                0.8,
                "--json",
                output_json,
            )
            payload = json.loads(output_json.read_text())
            sample = payload["samples"][0]
            self.assertEqual(sample["single_recalled_references"], 1)
            self.assertEqual(sample["union_recalled_references"], 1)
            self.assertEqual(sample["reciprocal_complete_references"], 1)
            self.assertAlmostEqual(sample["reference_aa_coverage_pct"], 90.0)


if __name__ == "__main__":
    unittest.main()
