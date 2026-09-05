#!/usr/bin/env python3
"""Stage 18: re-enter validated short rare fragments into the multi-k graph.

Root-cause hypothesis
---------------------
The production multi-k pipeline projects a stage to the next k using virtual
paired reads only for contigs >=500 bp. Stage10's validated rare rescue is
almost entirely 200--499 bp, so those fragments are appended after assembly but
never get a chance to recruit raw singleton k-mers and grow at k31/k41/k55.

This experiment changes only that handoff mechanism. Stage10 strict additions
are already cross-k validated; they are reintroduced as *single-copy* virtual
evidence. Virtual source intervals are non-overlapping, and a virtual pair is
rejected if any target-k k-mer is repeated inside that pair or has appeared in
another accepted virtual pair. Therefore the synthetic library alone cannot
satisfy --min-count=2 for a target-k k-mer; it can only promote sequence that
also has non-synthetic support.

After each target-k assembly, only contigs carrying a locus-unique Stage10 seed
signature are allowed to seed the next k. Final sequence is add-only on top of
Stage10 and must be seed-connected, Stage10-novel, and independently present in
at least two target-k assemblies.

No reference is used.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import resource
import subprocess
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import low_abundance_rescue as lr
import stage14_amplified_methods as s14
import stage16_root_cause as s16

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def run(cmd: list[object], *, env: dict[str, str] | None = None) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True, env=env)
    return time.monotonic() - started


@dataclass
class ConnectedRecord:
    name: str
    seq: str
    seed_id: int
    seed_hits: int
    second_hits: int


def virtual_pair_sequences(
    seq: str,
    *,
    read_length: int = 91,
    desired_insert: int = 220,
) -> Iterable[tuple[int, str, str]]:
    """Yield disjoint source intervals as paired virtual reads."""
    if desired_insert < 2 * read_length:
        raise ValueError("desired_insert must be >= 2*read_length")
    pos = 0
    while len(seq) - pos >= 2 * read_length:
        insert = min(desired_insert, len(seq) - pos)
        if insert < 2 * read_length:
            break
        left = seq[pos : pos + read_length]
        right_start = pos + insert - read_length
        right_forward = seq[right_start : right_start + read_length]
        if len(left) != read_length or len(right_forward) != read_length:
            break
        yield pos, left, rc(right_forward)
        pos += insert


def target_kmer_counts_from_pair(left: str, right_rc: str, k: int) -> Counter[str]:
    result: Counter[str] = Counter(lr.kmers(left, k))
    result.update(lr.kmers(rc(right_rc), k))
    return result


def write_single_copy_virtual_pairs(
    source_fasta: Path,
    out1: Path,
    out2: Path,
    *,
    target_k: int,
    read_length: int = 91,
    desired_insert: int = 220,
) -> dict[str, int | float]:
    out1.parent.mkdir(parents=True, exist_ok=True)
    seen_sequences: set[str] = set()
    used_target_kmers: set[str] = set()
    source_records = 0
    source_bases = 0
    accepted_pairs = 0
    accepted_read_bases = 0
    too_short_records = 0
    duplicate_pair_rejects = 0
    internal_repeat_rejects = 0
    n_pair_rejects = 0
    with gzip.open(out1, "wt", compresslevel=3) as left_out, gzip.open(
        out2, "wt", compresslevel=3
    ) as right_out:
        for _name, seq0 in lr.fasta_records(source_fasta):
            seq = seq0.upper()
            canonical = lr.canonical(seq)
            if canonical in seen_sequences:
                continue
            seen_sequences.add(canonical)
            source_records += 1
            source_bases += len(seq)
            emitted = 0
            for pos, left, right in virtual_pair_sequences(
                seq,
                read_length=read_length,
                desired_insert=desired_insert,
            ):
                if "N" in left or "N" in right:
                    n_pair_rejects += 1
                    continue
                counts = target_kmer_counts_from_pair(left, right, target_k)
                if not counts:
                    continue
                if max(counts.values()) > 1:
                    internal_repeat_rejects += 1
                    continue
                kmers = set(counts)
                if kmers & used_target_kmers:
                    duplicate_pair_rejects += 1
                    continue
                used_target_kmers.update(kmers)
                ident = f"handoff_k{target_k}_{source_records:06d}_{pos:05d}"
                qual = "I" * read_length
                left_out.write(f"@{ident}/1\n{left}\n+\n{qual}\n")
                right_out.write(f"@{ident}/2\n{right}\n+\n{qual}\n")
                accepted_pairs += 1
                accepted_read_bases += 2 * read_length
                emitted += 1
            if emitted == 0:
                too_short_records += int(len(seq) < 2 * read_length)
    return {
        "source_records": source_records,
        "source_bases": source_bases,
        "accepted_virtual_pairs": accepted_pairs,
        "accepted_virtual_read_bases": accepted_read_bases,
        "unique_virtual_target_kmers": len(used_target_kmers),
        "duplicate_pair_rejects": duplicate_pair_rejects,
        "internal_repeat_rejects": internal_repeat_rejects,
        "n_pair_rejects": n_pair_rejects,
        "too_short_records": too_short_records,
        "virtual_to_source_base_ratio": accepted_read_bases / max(1, source_bases),
    }


def concat_gzip(inputs: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as out:
        for path in inputs:
            with path.open("rb") as handle:
                while True:
                    chunk = handle.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)


def assemble_target_k(
    bridgeasm: Path,
    read1: Path,
    read2: Path,
    outdir: Path,
    k: int,
    threads: int,
) -> float:
    mercy = {31: 16, 41: 12, 55: 8}[k]
    env = os.environ.copy()
    env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION", None)
    env.pop("BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY", None)
    env.pop("BRIDGEASM_MATE_TERMINAL_MERCY_KMERS", None)
    return run(
        [
            bridgeasm,
            "assemble",
            "-1",
            read1,
            "-2",
            read2,
            "-o",
            outdir,
            "-k",
            k,
            "--min-count",
            2,
            "--mercy-max-kmers",
            mercy,
            "--mercy-min-support",
            1,
            "--mercy-min-quality",
            25,
            "--min-read-support",
            2,
            "--min-pair-support",
            2,
            "--min-primary-support",
            5,
            "--primary-dominance",
            0.75,
            "--threaded-path-cover",
            "--major-path-cover",
            "--path-cover-secondary-dominance",
            0.25,
            "--min-contig-length",
            200,
            "--threads",
            threads,
        ],
        env=env,
    )


def assign_sequence_to_seed(
    seq: str,
    signatures: dict[str, int],
    *,
    k: int = 21,
    min_hits: int = 4,
    margin: int = 2,
) -> tuple[int | None, int, int]:
    scores: Counter[int] = Counter()
    for mer in set(lr.kmers(seq, k)):
        sid = signatures.get(mer)
        if sid is not None:
            scores[sid] += 1
    ranked = scores.most_common(2)
    if not ranked:
        return None, 0, 0
    best_sid, best_hits = ranked[0]
    second_hits = ranked[1][1] if len(ranked) > 1 else 0
    if best_hits < min_hits or best_hits < second_hits + margin:
        return None, best_hits, second_hits
    return best_sid, best_hits, second_hits


def filter_seed_connected(
    inputs: list[Path],
    signatures: dict[str, int],
    output: Path,
) -> tuple[list[ConnectedRecord], dict[str, int]]:
    records: list[ConnectedRecord] = []
    seen: set[str] = set()
    input_records = 0
    ambiguous_or_weak = 0
    for path in inputs:
        if not path.exists() or path.stat().st_size == 0:
            continue
        for name, seq0 in lr.fasta_records(path):
            input_records += 1
            seq = lr.canonical(seq0)
            if seq in seen:
                continue
            seen.add(seq)
            sid, hits, second = assign_sequence_to_seed(seq, signatures)
            if sid is None:
                ambiguous_or_weak += 1
                continue
            records.append(ConnectedRecord(name, seq, sid, hits, second))
    records.sort(key=lambda item: (item.seed_id, -item.seed_hits, -len(item.seq), item.seq))
    s14.write_fasta(
        (
            (f"seed{item.seed_id}_hits{item.seed_hits}_{idx:06d}", item.seq)
            for idx, item in enumerate(records, 1)
        ),
        output,
    )
    return records, {
        "input_records": input_records,
        "unique_sequences": len(seen),
        "seed_connected_records": len(records),
        "seed_connected_bases": sum(len(item.seq) for item in records),
        "connected_seed_loci": len({item.seed_id for item in records}),
        "ambiguous_or_weak_records": ambiguous_or_weak,
    }


def profile(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.exists() else {}


def build_short_handoff(
    scripts: Path,
    bridgeasm: Path,
    pipeline_dir: Path,
    strict_baseline: Path,
    backbone: Path,
    read1: Path,
    read2: Path,
    threads: int,
    timings: dict[str, float],
) -> tuple[Path, Path, dict[str, object]]:
    stage10 = pipeline_dir / "stage10_multik_rescue"
    seed_fasta = stage10 / "multik_strict_additions.fasta"
    root = pipeline_dir / "stage18_short_handoff"
    root.mkdir(parents=True, exist_ok=True)

    seeds, signatures, _sets, signature_stats = s16.build_all_seed_signatures(
        seed_fasta, backbone, k=21, min_signature_kmers=4
    )
    source = seed_fasta
    assembly_inputs: dict[int, list[Path]] = {}
    stage_stats: dict[str, object] = {}

    for k in (31, 41, 55):
        stage = root / f"k{k}"
        stage.mkdir(parents=True, exist_ok=True)
        virtual1 = stage / "handoff_R1.fastq.gz"
        virtual2 = stage / "handoff_R2.fastq.gz"
        vstats = write_single_copy_virtual_pairs(
            source,
            virtual1,
            virtual2,
            target_k=k,
            read_length=91,
            desired_insert=220,
        )
        aug1 = stage / "aug_R1.fastq.gz"
        aug2 = stage / "aug_R2.fastq.gz"
        concat_gzip([read1, virtual1], aug1)
        concat_gzip([read2, virtual2], aug2)
        asm = stage / "assembly"
        timings[f"short_handoff_k{k}"] = assemble_target_k(
            bridgeasm, aug1, aug2, asm, k, threads
        )
        connected = stage / "seed_connected.fasta"
        _records, cstats = filter_seed_connected(
            [asm / "primary_contigs.fasta", asm / "haplotigs.fasta"],
            signatures,
            connected,
        )
        assembly_inputs[k] = [connected]
        stage_stats[f"k{k}"] = {
            "virtual": vstats,
            "assembly_profile": profile(asm / "run_profile.json"),
            "connected": cstats,
        }
        source = connected

    evidence = s16.local_evidence(
        assembly_inputs,
        seeds,
        strict_baseline,
        Counter(),
    )
    selected = s16.select_local_evidence(evidence, strict_baseline)
    s16.write_local_evidence(evidence, root / "short_handoff_evidence.tsv")
    s16.write_local_evidence(selected, root / "short_handoff_selected.tsv")
    additions = root / "short_handoff_additions.fasta"
    s14.write_fasta(
        (
            (f"stage18_seed{item.seed_id}_k{item.k}_{idx:06d}", item.seq)
            for idx, item in enumerate(selected, 1)
        ),
        additions,
    )
    final = s14.make_bridge_candidate(
        scripts,
        strict_baseline,
        [additions],
        root / "candidate_short_handoff",
        timings,
        min_overlap=81,
    )
    return final, additions, {
        **signature_stats,
        "stages": stage_stats,
        "evidence_records": len(evidence),
        "selected_records": len(selected),
        "selected_bases": sum(len(item.seq) for item in selected),
        "selected_seed_loci": len({item.seed_id for item in selected}),
        "selected_fresh31": sum(item.fresh31 for item in selected),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    scripts = Path(__file__).resolve().parent
    pipeline = args.pipeline_dir
    stage10 = pipeline / "stage10_multik_rescue"
    strict_baseline = stage10 / "candidate_multik_strict" / "primary_contigs.fasta"
    backbone = pipeline / "bridge_backbone.fasta"
    required = [
        args.bridgeasm,
        args.read1,
        args.read2,
        strict_baseline,
        backbone,
        stage10 / "multik_strict_additions.fasta",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing Stage18 inputs: " + ", ".join(missing))

    started = time.monotonic()
    timings: dict[str, float] = {}
    final, additions, method_stats = build_short_handoff(
        scripts,
        args.bridgeasm,
        pipeline,
        strict_baseline,
        backbone,
        args.read1,
        args.read2,
        args.threads,
        timings,
    )
    root = pipeline / "stage18_short_handoff"
    stats = {
        "pipeline": "bridge-stage18-short-handoff-v1",
        "baseline": str(strict_baseline),
        "policy": {
            "reference_free": True,
            "metric_targets": False,
            "source": "Stage10 strict cross-k validated rare additions",
            "handoff": "single-copy non-overlapping virtual evidence k31->k41->k55",
            "synthetic_solid_kmer_invariant": "target-k synthetic multiplicity <=1 globally",
            "next_stage_filter": "locus-unique Stage10 seed signatures",
            "final_validation": "seed-connected + Stage10-novel + cross-target-k recurrence",
            "sequence_join": "exact overlap >=81 bp",
        },
        "method": method_stats,
        "outputs": {
            "final": str(final),
            "additions": str(additions),
        },
        "timings_seconds": timings,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (root / "stage18_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
