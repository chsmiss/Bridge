#!/usr/bin/env python3
from __future__ import annotations

import gzip
import random
import tempfile
from collections import Counter
from pathlib import Path

import low_abundance_rescue as lr
import stage16_root_cause as s16
import stage18_short_handoff as s18


def random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def fastq_seqs(path: Path) -> list[str]:
    out: list[str] = []
    with gzip.open(path, "rt") as handle:
        while True:
            h = handle.readline()
            if not h:
                break
            out.append(handle.readline().strip())
            handle.readline()
            handle.readline()
    return out


def test_virtual_source_intervals_are_disjoint() -> None:
    seq = random_dna(660, 1)
    pairs = list(s18.virtual_pair_sequences(seq, read_length=91, desired_insert=220))
    assert [pos for pos, _l, _r in pairs] == [0, 220, 440]
    assert len(pairs) == 3


def test_virtual_target_kmers_have_global_copy_one() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "source.fasta"
        source.write_text(">a\n" + random_dna(440, 2) + "\n")
        r1 = root / "r1.fastq.gz"
        r2 = root / "r2.fastq.gz"
        stats = s18.write_single_copy_virtual_pairs(
            source, r1, r2, target_k=31, read_length=91, desired_insert=220
        )
        assert stats["accepted_virtual_pairs"] == 2
        counts: Counter[str] = Counter()
        for seq in fastq_seqs(r1):
            counts.update(lr.kmers(seq, 31))
        for seq in fastq_seqs(r2):
            counts.update(lr.kmers(s18.rc(seq), 31))
        assert counts
        assert max(counts.values()) == 1


def test_repeated_virtual_kmer_pair_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "source.fasta"
        source.write_text(">repeat\n" + ("ACGT" * 55) + "\n")
        r1 = root / "r1.fastq.gz"
        r2 = root / "r2.fastq.gz"
        stats = s18.write_single_copy_virtual_pairs(
            source, r1, r2, target_k=31, read_length=91, desired_insert=220
        )
        assert stats["accepted_virtual_pairs"] == 0
        assert stats["internal_repeat_rejects"] >= 1


def test_locus_unique_signature_assignment() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        backbone = root / "backbone.fasta"
        seeds = root / "seeds.fasta"
        backbone.write_text(">b\n" + random_dna(400, 10) + "\n")
        seq0 = random_dna(260, 20)
        seq1 = random_dna(260, 30)
        seeds.write_text(f">s0\n{seq0}\n>s1\n{seq1}\n")
        _records, signatures, _sets, stats = s16.build_all_seed_signatures(
            seeds, backbone, k=21, min_signature_kmers=4
        )
        assert stats["signature_eligible_seeds"] == 2
        sid, hits, second = s18.assign_sequence_to_seed(seq0 + random_dna(80, 40), signatures)
        assert sid == 0
        assert hits >= 4
        assert second == 0


def main() -> None:
    tests = [
        test_virtual_source_intervals_are_disjoint,
        test_virtual_target_kmers_have_global_copy_one,
        test_repeated_virtual_kmer_pair_is_rejected,
        test_locus_unique_signature_assignment,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)


if __name__ == "__main__":
    main()
