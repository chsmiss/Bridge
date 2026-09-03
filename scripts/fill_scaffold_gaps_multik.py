#!/usr/bin/env python3
"""Fill pair-scaffold gaps only when multiple local-k assemblies agree.

This is deliberately stricter than ordinary scaffolding: pair links may propose
positive gaps, but an N-gap only becomes sequence if at least two independent
local de Bruijn k values recover exactly the same fill. Any unresolved N-gap is
split before emission, so reported contig/NA50 gains are sequence-resolved and
not caused by artificial N scaffolds.
"""
from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

import fill_scaffold_gaps as fg


def consensus_fill(fills: dict[int, str | None], min_consensus: int) -> tuple[str | None, list[int]]:
    groups: defaultdict[str, list[int]] = defaultdict(list)
    for k, seq in fills.items():
        if seq is not None:
            groups[seq].append(k)
    if not groups:
        return None, []
    ranked = sorted(groups.items(), key=lambda item: (-len(item[1]), -len(item[0]), item[0]))
    seq, ks = ranked[0]
    if len(ks) < min_consensus:
        return None, []
    return seq, sorted(ks)


def split_unresolved(records: list[tuple[str, str]], min_length: int) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for name, seq in records:
        parts = [part for part in re.split(r"N+", seq.upper()) if len(part) >= min_length]
        for index, part in enumerate(parts, 1):
            out.append((f"{name}_part{index}", part))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("scaffolds", type=Path)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--anchor-k", type=int, default=31)
    ap.add_argument("--local-ks", default="17,21,25")
    ap.add_argument("--min-consensus", type=int, default=2)
    ap.add_argument("--flank", type=int, default=220)
    ap.add_argument("--dominance", type=float, default=0.68)
    ap.add_argument("--max-reads-per-gap", type=int, default=1500)
    ap.add_argument("--min-length", type=int, default=200)
    args = ap.parse_args()

    ks = sorted({int(x) for x in args.local_ks.split(",") if x.strip()})
    if len(ks) < args.min_consensus or any(k < 13 or k > args.anchor_k for k in ks):
        raise SystemExit("invalid --local-ks/--min-consensus")

    recs = fg.fasta(args.scaffolds)
    gaps: list[dict[str, object]] = []
    anchor_index: defaultdict[str, list[int]] = defaultdict(list)
    for rid, (_name, seq) in enumerate(recs):
        for match in re.finditer(r"N+", seq):
            ll = max(0, match.start() - args.flank)
            rr = min(len(seq), match.end() + args.flank)
            left = seq[ll : match.start()]
            right = seq[match.end() : rr]
            if len(left) < args.anchor_k or len(right) < args.anchor_k:
                continue
            gid = len(gaps)
            gaps.append(
                {
                    "rid": rid,
                    "start": match.start(),
                    "end": match.end(),
                    "left": left,
                    "right": right,
                    "reads": [],
                }
            )
            for q in fg.anchor_kmers(left[-min(len(left), 120) :], args.anchor_k):
                anchor_index[q].append(gid)
            for q in fg.anchor_kmers(right[: min(len(right), 120)], args.anchor_k):
                anchor_index[q].append(gid)

    for (_id1, s1), (_id2, s2) in fg.fastq_pairs(args.read1, args.read2):
        hit: set[int] = set()
        for seq in (s1, s2):
            for pos in range(len(seq) - args.anchor_k + 1):
                q = fg.canon(seq[pos : pos + args.anchor_k])
                hit.update(anchor_index.get(q, ()))
        for gid in hit:
            buf = gaps[gid]["reads"]
            assert isinstance(buf, list)
            if len(buf) < args.max_reads_per_gap:
                buf.extend((s1, s2))

    replacements: defaultdict[int, list[tuple[int, int, str]]] = defaultdict(list)
    rows: list[tuple[object, ...]] = []
    for gid, gap in enumerate(gaps):
        left = str(gap["left"])
        right = str(gap["right"])
        reads = gap["reads"]
        assert isinstance(reads, list)
        expected = int(gap["end"]) - int(gap["start"])
        attempts: dict[int, str | None] = {}
        for k in ks:
            attempts[k] = fg.local_path(
                left,
                right,
                reads,
                k,
                expected,
                expected + 600,
                args.dominance,
            )
        fill, support_ks = consensus_fill(attempts, args.min_consensus)
        status = "consensus_filled" if fill is not None else "unresolved"
        if fill is not None:
            replacements[int(gap["rid"])].append(
                (int(gap["start"]), int(gap["end"]), fill)
            )
        rows.append(
            (
                gid,
                recs[int(gap["rid"])][0],
                int(gap["start"]),
                int(gap["end"]),
                expected,
                len(reads),
                status,
                0 if fill is None else len(fill),
                ",".join(map(str, support_ks)),
                ",".join(f"k{k}:{0 if seq is None else len(seq)}" for k, seq in sorted(attempts.items())),
            )
        )

    replaced: list[tuple[str, str]] = []
    for rid, (name, seq0) in enumerate(recs):
        seq = seq0
        for start, end, fill in sorted(replacements[rid], reverse=True):
            seq = seq[:start] + fill + seq[end:]
        replaced.append((name, seq))
    final = split_unresolved(replaced, args.min_length)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for index, (_name, seq) in enumerate(final, 1):
            handle.write(f">multik_gap_refined_{index:07d} len={len(seq)}\n")
            for pos in range(0, len(seq), 80):
                handle.write(seq[pos : pos + 80] + "\n")
    with args.report.open("w") as handle:
        handle.write(
            "gap_id\tscaffold\tstart\tend\testimated_gap\tlocal_reads\tstatus"
            "\tfilled_bases\tconsensus_ks\tattempts\n"
        )
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")

    status_counts = Counter(row[6] for row in rows)
    print(f"gaps\t{len(gaps)}")
    print(f"consensus_filled\t{status_counts.get('consensus_filled', 0)}")
    print(f"unresolved\t{status_counts.get('unresolved', 0)}")
    print(f"output_records\t{len(final)}")


if __name__ == "__main__":
    main()
