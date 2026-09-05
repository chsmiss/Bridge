#!/usr/bin/env python3
from __future__ import annotations

import random
import subprocess
import tempfile
from pathlib import Path


def seq(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def fasta_sequences(path: Path) -> list[str]:
    out: list[str] = []
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if chunks:
                out.append("".join(chunks))
                chunks = []
        else:
            chunks.append(raw.strip())
    if chunks:
        out.append("".join(chunks))
    return out


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        backbone = seq(1, 600)
        novel = seq(2, 55)
        short_tail = backbone[-64:] + novel
        inp = td / "in.fasta"
        inp.write_text(f">backbone\n{backbone}\n>tail\n{short_tail}\n")
        script = Path(__file__).resolve().parent / "postprocess_segment_grafts.py"
        subprocess.run(
            [
                "python3",
                str(script),
                str(inp),
                "--output-dir",
                str(td / "out"),
                "--short-min-length",
                "64",
                "--final-min-length",
                "200",
                "--min-overlap",
                "31",
            ],
            check=True,
        )
        outputs = fasta_sequences(td / "out" / "primary_contigs.fasta")
        expected = backbone + novel
        rc = str.maketrans("ACGT", "TGCA")
        expected_rc = expected.translate(rc)[::-1]
        assert expected in outputs or expected_rc in outputs, outputs
        assert all(len(item) >= 200 for item in outputs)
    print("short segment graft synthetic test: ok")


if __name__ == "__main__":
    main()
