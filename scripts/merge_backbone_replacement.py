#!/usr/bin/env python3
"""Reference-free replacement merge for a new graph backbone and recovery contigs.

The graph backbone is authoritative for sequence regions it covers. Recovery
contigs whose canonical k-mers are mostly represented by the backbone are
removed before union, preventing partial-overlap duplicates from inflating Dup.
Recovery contigs with genuinely novel sequence are retained.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

BASE = {65: 0, 67: 1, 71: 2, 84: 3, 97: 0, 99: 1, 103: 2, 116: 3}


def fasta(path: Path):
    name = None
    chunks = []
    with path.open() as handle:
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
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def canonical_keys(seq: str, k: int):
    mask = (1 << (2 * k)) - 1
    forward = reverse = valid = 0
    shift = 2 * (k - 1)
    for byte in seq.encode("ascii", "ignore"):
        value = BASE.get(byte)
        if value is None:
            forward = reverse = valid = 0
            continue
        forward = ((forward << 2) | value) & mask
        reverse = (reverse >> 2) | ((3 - value) << shift)
        valid += 1
        if valid >= k:
            yield min(forward, reverse)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", type=Path, required=True)
    ap.add_argument("--recovery", type=Path, required=True)
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("--replace-fraction", type=float, default=0.85)
    ap.add_argument("--min-informative-kmers", type=int, default=20)
    args = ap.parse_args()

    backbone = list(fasta(args.backbone))
    recovery = list(fasta(args.recovery))
    backbone_keys = set()
    for _, seq in backbone:
        backbone_keys.update(canonical_keys(seq, args.k))

    kept = []
    rows = []
    removed_bases = 0
    for name, seq in recovery:
        keys = list(canonical_keys(seq, args.k))
        informative = len(keys)
        represented = sum(key in backbone_keys for key in keys)
        fraction = represented / max(1, informative)
        replace = (
            informative >= args.min_informative_kmers
            and fraction >= args.replace_fraction
        )
        if replace:
            removed_bases += len(seq)
        else:
            kept.append((name, seq))
        rows.append((name, len(seq), informative, represented, fraction, int(replace)))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        index = 0
        for source, records in (("backbone", backbone), ("recovery", kept)):
            for name, seq in records:
                index += 1
                out.write(
                    f">{source}_{index:07d} original={name} len={len(seq)}\n"
                )
                for start in range(0, len(seq), 80):
                    out.write(seq[start : start + 80] + "\n")

    with args.report.open("w") as out:
        out.write(
            "contig\tlength\tinformative_kmers\trepresented_kmers\t"
            "represented_fraction\treplaced\n"
        )
        for row in rows:
            out.write(
                "\t".join(map(str, row[:4]))
                + f"\t{row[4]:.6f}\t{row[5]}\n"
            )

    stats = {
        "backbone_records": len(backbone),
        "backbone_bases": sum(len(seq) for _, seq in backbone),
        "recovery_records": len(recovery),
        "recovery_bases": sum(len(seq) for _, seq in recovery),
        "recovery_replaced_records": len(recovery) - len(kept),
        "recovery_replaced_bases": removed_bases,
        "recovery_kept_records": len(kept),
        "recovery_kept_bases": sum(len(seq) for _, seq in kept),
        "replace_fraction": args.replace_fraction,
        "k": args.k,
        "output_records": len(backbone) + len(kept),
        "output_bases": sum(len(seq) for _, seq in backbone)
        + sum(len(seq) for _, seq in kept),
    }
    print(json.dumps(stats, sort_keys=True))


if __name__ == "__main__":
    main()
