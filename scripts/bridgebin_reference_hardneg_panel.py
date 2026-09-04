#!/usr/bin/env python3
"""Build a deterministic hard-negative BridgeBin training panel from reference genomes.

Input is a directory containing one subdirectory per biological group/species.  Each
subdirectory may contain one or more FASTA/FNA files.  The script samples long windows
from the references, writes a compact contig panel and truth table, and deliberately
assigns nearly identical synthetic multi-sample coverage profiles to every group.  This
prevents the pair head from solving the curriculum by abundance alone and forces DNA,
gene or protein identity to carry the difficult decision.

This panel is for training/validation curriculum only; it is not a benchmark score.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
from pathlib import Path
from typing import Iterable, Iterator, List, Sequence, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--fragment-bp", type=int, default=6000)
    p.add_argument("--fragments-per-group", type=int, default=8)
    p.add_argument("--samples", type=int, default=4)
    p.add_argument("--seed", type=int, default=43)
    p.add_argument("--coverage-noise", type=float, default=0.015)
    return p.parse_args()


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
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


def fasta_files(root: Path) -> List[Path]:
    suffixes = {".fa", ".fasta", ".fna", ".fas"}
    return sorted(path for path in root.rglob("*") if path.is_file() and path.suffix.lower() in suffixes)


def stable_seed(text: str, seed: int) -> int:
    digest = hashlib.blake2b(f"{seed}\0{text}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little")


def candidate_fragments(group_dir: Path, fragment_bp: int) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for fasta in fasta_files(group_dir):
        for record, sequence in read_fasta(fasta):
            if len(sequence) < fragment_bp:
                continue
            # Use non-overlapping windows to avoid trivially duplicated positives.
            for start in range(0, len(sequence) - fragment_bp + 1, fragment_bp):
                piece = sequence[start : start + fragment_bp]
                if piece.count("N") / max(1, len(piece)) <= 0.02:
                    out.append((record, start, piece))
    return out


def safe_name(value: str) -> str:
    return "_".join(piece for piece in value.replace("-", "_").split() if piece)


def main() -> int:
    args = parse_args()
    if args.fragment_bp < 1000 or args.fragments_per_group < 2 or args.samples < 1:
        raise SystemExit("fragment size/group count/samples are too small")
    if not 0.0 <= args.coverage_noise <= 0.2:
        raise SystemExit("--coverage-noise must be in [0,0.2]")

    groups = sorted(path for path in args.input_dir.iterdir() if path.is_dir())
    if len(groups) < 3:
        raise SystemExit("need at least three group directories")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    selected: List[Tuple[str, str, str, int, str]] = []
    for group_dir in groups:
        group = safe_name(group_dir.name)
        candidates = candidate_fragments(group_dir, args.fragment_bp)
        if len(candidates) < args.fragments_per_group:
            raise SystemExit(
                f"{group_dir}: only {len(candidates)} usable fragments, need {args.fragments_per_group}"
            )
        rng = random.Random(stable_seed(group, args.seed))
        indices = sorted(rng.sample(range(len(candidates)), args.fragments_per_group))
        for serial, index in enumerate(indices, start=1):
            record, start, sequence = candidates[index]
            contig = f"{group}__ref{serial:03d}"
            selected.append((contig, group, record, start, sequence))

    with (args.output_dir / "contigs.fa").open("w", encoding="utf-8") as handle:
        for contig, _group, _record, _start, sequence in selected:
            handle.write(f">{contig}\n")
            for pos in range(0, len(sequence), 80):
                handle.write(sequence[pos : pos + 80] + "\n")

    with (args.output_dir / "truth.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "genome", "length", "eligible", "record", "start", "end"])
        for contig, group, record, start, sequence in selected:
            writer.writerow([contig, group, len(sequence), 1, record, start, start + len(sequence)])

    # Every group gets the same abundance backbone.  Tiny deterministic contig-level
    # noise prevents exact floating-point duplicates while leaving coverage effectively
    # useless for distinguishing biological groups.
    base = [30.0, 12.0, 25.0, 7.0, 18.0, 10.0, 22.0, 14.0][: args.samples]
    if args.samples > len(base):
        base.extend(15.0 + 2.0 * i for i in range(args.samples - len(base)))
    sample_names = [f"sample{i+1}" for i in range(args.samples)]
    with (args.output_dir / "coverage.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", *sample_names])
        for contig, group, _record, _start, _sequence in selected:
            rng = random.Random(stable_seed(contig, args.seed))
            values = [depth * (1.0 + rng.uniform(-args.coverage_noise, args.coverage_noise)) for depth in base]
            writer.writerow([contig, *[f"{value:.6f}" for value in values]])

    print(
        f"bridgebin-hardneg-panel: groups={len(groups)} contigs={len(selected)} "
        f"fragment_bp={args.fragment_bp} samples={args.samples}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
