#!/usr/bin/env python3
from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import low_abundance_rescue as lr
import stage19_trusted_fragment_handoff as s19


def read_fastq(path: Path) -> list[tuple[str, str, str]]:
    rows=[]
    with gzip.open(path,'rt') as h:
        while True:
            name=h.readline().strip()
            if not name: break
            seq=h.readline().strip(); h.readline(); qual=h.readline().strip()
            rows.append((name,seq,qual))
    return rows


def test_all_priors_are_one_physical_fragment() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        fa=root/'p.fa'
        a='ACGTGCAATGCTAGCTACGATCGTACGTTAGC'*3
        b='TTGACCGTATCGGATCCGATGCTAGCATTCGA'*3
        fa.write_text(f'>a\n{a}\n>b\n{b}\n')
        r1=root/'r1.fq.gz'; r2=root/'r2.fq.gz'
        stats=s19.write_trusted_fragment(fa,r1,r2,target_k=31,synthetic_phred=20)
        left=read_fastq(r1); right=read_fastq(r2)
        assert len(left)==1 and len(right)==1
        assert stats['synthetic_fragments']==1
        assert ('N'*31) in left[0][1]
        assert right[0][1]=='N'
        assert set(left[0][2])=={'5'}  # ASCII 53 = Phred 20


def test_separator_prevents_cross_seed_kmers() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        fa=root/'p.fa'
        a='A'*40+'C'*40
        b='G'*40+'T'*40
        fa.write_text(f'>a\n{a}\n>b\n{b}\n')
        r1=root/'r1.fq.gz'; r2=root/'r2.fq.gz'
        s19.write_trusted_fragment(fa,r1,r2,target_k=31)
        joined=read_fastq(r1)[0][1]
        observed=set(lr.kmers(joined,31))
        expected=set(lr.kmers(a,31)) | set(lr.kmers(b,31))
        assert observed==expected


def test_duplicate_prior_kmer_still_has_one_fragment_semantically() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        fa=root/'p.fa'
        motif='ACGTTGCAACGTAGCTTGCAACGTTGCAACGTAGCTTGCA'
        fa.write_text(f'>a\n{motif}\n>b\n{motif}AA\n')
        r1=root/'r1.fq.gz'; r2=root/'r2.fq.gz'
        stats=s19.write_trusted_fragment(fa,r1,r2,target_k=31)
        # There is one FASTQ pair total. BridgeAsm deduplicates fragment_keys
        # before incrementing fragment_count, so repeated observations remain
        # one physical-fragment support event.
        assert len(read_fastq(r1))==1
        assert stats['max_prior_fragment_support']==1
        assert stats['target_kmer_observations'] >= stats['distinct_target_kmers']


def main() -> None:
    tests=[
        test_all_priors_are_one_physical_fragment,
        test_separator_prevents_cross_seed_kmers,
        test_duplicate_prior_kmer_still_has_one_fragment_semantically,
    ]
    for test in tests:
        test(); print('PASS',test.__name__)


if __name__=='__main__': main()
