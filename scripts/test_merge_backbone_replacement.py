#!/usr/bin/env python3
"""Synthetic checks for segment-level backbone replacement."""
from __future__ import annotations

import random
import subprocess
import tempfile
from pathlib import Path


def random_seq(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def read_fasta(path: Path) -> list[str]:
    records: list[str] = []
    name = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if name is not None:
                records.append("".join(chunks))
            name = raw[1:]
            chunks = []
        else:
            chunks.append(raw.strip())
    if name is not None:
        records.append("".join(chunks))
    return records


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        backbone = random_seq(1, 1200)
        novel_tail = random_seq(2, 120)
        novel_mid = random_seq(3, 45)
        mostly_novel = random_seq(4, 350)
        recovery_tail = backbone[:1000] + novel_tail
        recovery_exact = backbone[150:950]
        recovery_mid = backbone[100:600] + novel_mid + backbone[645:1050]

        (out / "backbone.fasta").write_text(f">backbone\n{backbone}\n")
        (out / "recovery.fasta").write_text(
            f">tail\n{recovery_tail}\n"
            f">exact\n{recovery_exact}\n"
            f">middle\n{recovery_mid}\n"
            f">novel\n{mostly_novel}\n"
        )
        script = Path(__file__).resolve().parent / "merge_backbone_replacement.py"
        subprocess.run(
            [
                "python3",
                str(script),
                "--backbone",
                str(out / "backbone.fasta"),
                "--recovery",
                str(out / "recovery.fasta"),
                "-o",
                str(out / "merged.fasta"),
                "--report",
                str(out / "report.tsv"),
                "-k",
                "21",
                "--replace-fraction",
                "0.80",
                "--segment-anchor-bases",
                "50",
                "--min-novel-kmers",
                "4",
                "--min-segment-length",
                "120",
            ],
            check=True,
        )
        seqs = read_fasta(out / "merged.fasta")
        assert backbone in seqs
        assert recovery_exact not in seqs
        assert mostly_novel in seqs
        assert recovery_tail not in seqs
        assert recovery_mid not in seqs
        assert any(
            seq.endswith(novel_tail) and len(seq) < len(recovery_tail) for seq in seqs
        )
        assert any(novel_mid in seq and len(seq) < len(recovery_mid) for seq in seqs)
        report = (out / "report.tsv").read_text()
        assert "keep_segment" in report
        assert "replace_full" in report
    print("segment-level replacement synthetic test: ok")


if __name__ == "__main__":
    main()
