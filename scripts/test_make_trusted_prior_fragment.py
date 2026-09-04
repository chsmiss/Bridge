#!/usr/bin/env python3
from __future__ import annotations

import gzip
import random
import tempfile
from pathlib import Path

import low_abundance_rescue as lr
import make_trusted_prior_fragment as tp


def dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def fastq_records(path: Path) -> list[tuple[str, str, str]]:
    rows = []
    with gzip.open(path, "rt") as handle:
        while True:
            name = handle.readline().strip()
            if not name:
                break
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            rows.append((name, seq, qual))
    return rows


def test_short_validated_prior_is_not_cut_at_500() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "prior.fa"
        source.write_text(">short\n" + dna(220, 1) + "\n")
        seqs, stats = tp.build_prior_sequences([source], min_length=200)
        assert len(seqs) == 1
        assert len(seqs[0]) == 220
        assert stats["short_records"] == 0


def test_reverse_complement_duplicate_counts_once() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        seq = dna(260, 2)
        rc = lr.canonical(seq)
        # canonical(seq) is either seq or its RC; writing both orientations
        # must still leave one structural prior record.
        comp = str.maketrans("ACGT", "TGCA")
        reverse = seq.translate(comp)[::-1]
        source = root / "prior.fa"
        source.write_text(f">a\n{seq}\n>b\n{reverse}\n")
        seqs, stats = tp.build_prior_sequences([source], min_length=200)
        assert len(seqs) == 1
        assert stats["duplicate_records"] == 1
        assert lr.canonical(seqs[0]) == rc


def test_all_records_are_one_physical_pair_and_n_separated() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        seqs = [dna(230, 3), dna(250, 4)]
        r1 = root / "r1.fastq.gz"
        r2 = root / "r2.fastq.gz"
        stats = tp.write_one_fragment(seqs, r1, r2, target_k=31, phred=20)
        left = fastq_records(r1)
        right = fastq_records(r2)
        assert len(left) == 1 and len(right) == 1
        assert ("N" * 31) in left[0][1]
        assert right[0][1] == "N"
        assert stats["synthetic_physical_fragments"] == 1
        assert stats["max_synthetic_fragment_support_per_kmer"] == 1
        assert set(left[0][2]) == {"5"}


def test_separator_creates_no_cross_record_target_kmer() -> None:
    seqs = [dna(220, 5), dna(220, 6)]
    joined = ("N" * 31).join(seqs)
    observed = set(lr.kmers(joined, 31))
    expected = set(lr.kmers(seqs[0], 31)) | set(lr.kmers(seqs[1], 31))
    assert observed == expected


def main() -> None:
    tests = [
        test_short_validated_prior_is_not_cut_at_500,
        test_reverse_complement_duplicate_counts_once,
        test_all_records_are_one_physical_pair_and_n_separated,
        test_separator_creates_no_cross_record_target_kmer,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)


if __name__ == "__main__":
    main()
