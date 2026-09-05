#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Dict, Optional, Sequence

CODONS: Dict[str, str] = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def random_protein(rng: random.Random, length: int) -> str:
    return "M" + "".join(rng.choice(AA_ALPHABET) for _ in range(length - 1))


def back_translate(protein: str) -> str:
    return "".join(CODONS[residue] for residue in protein)


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def write_fasta(path: Path, records) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{wrap(sequence)}\n")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=43)
    parser.add_argument("--protein-aa", type=int, default=320)
    parser.add_argument("--read-nt", type=int, default=300)
    parser.add_argument("--step-nt", type=int, default=24)
    parser.add_argument("--copies", type=int, default=8)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    protein = random_protein(rng, args.protein_aa)
    dna = back_translate(protein)

    split = len(dna) // 2
    overlap = 31
    source = dna[:split]
    target = dna[split - overlap :]
    decoy_protein = random_protein(rng, (len(dna) - split + 2) // 3 + 20)
    decoy_suffix = back_translate(decoy_protein)[: len(dna) - split]
    decoy = dna[split - overlap : split] + decoy_suffix

    graph = args.output / "graph.gfa"
    graph.write_text(
        "\n".join(
            [
                "H\tVN:Z:1.0",
                f"S\tA\t{source}\tKC:f:40.0",
                f"S\tB\t{target}\tKC:f:40.0",
                f"S\tC\t{decoy}\tKC:f:40.0",
                f"L\tA\t+\tB\t+\t{overlap}M\tDR:i:2\tGR:i:0\tPE:i:2",
                f"L\tA\t+\tC\t+\t{overlap}M\tDR:i:2\tGR:i:0\tPE:i:2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    write_fasta(args.output / "reference_protein.faa", [("reference_protein", protein)])
    write_fasta(args.output / "expected_dna.fasta", [("expected_true_path", dna)])

    reads = []
    read_index = 0
    starts = list(range(0, max(1, len(dna) - args.read_nt + 1), args.step_nt))
    if not starts or starts[-1] != len(dna) - args.read_nt:
        starts.append(max(0, len(dna) - args.read_nt))
    for copy in range(args.copies):
        for start in starts:
            sequence = dna[start : start + args.read_nt]
            if len(sequence) < 90:
                continue
            read_index += 1
            if (copy + start // max(1, args.step_nt)) % 2:
                sequence = reverse_complement(sequence)
                orientation = "rc"
            else:
                orientation = "fw"
            reads.append((f"read_{read_index:06}_{orientation}_start_{start}", sequence))
    write_fasta(args.output / "reads.fasta", reads)

    manifest = {
        "seed": args.seed,
        "protein_aa": len(protein),
        "dna_nt": len(dna),
        "read_count": len(reads),
        "read_nt": args.read_nt,
        "step_nt": args.step_nt,
        "copies": args.copies,
        "split": split,
        "overlap": overlap,
        "true_edge": ["A", "B"],
        "decoy_edge": ["A", "C"],
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
