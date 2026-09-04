#!/usr/bin/env python3
from __future__ import annotations

import random

import graph_path_phaser as gp
import stage20_anchor_gap_thread as s20


def dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(length))


def bubble_graph() -> gp.Graph:
    g=gp.Graph(); g.k=31
    g.seqs={'A':dna(70,1),'B':dna(70,2),'C':dna(70,3),'D':dna(70,4)}
    g.coverage={u:1.0 for u in g.seqs}
    for u in g.seqs:
        g.out[u]=[]; g.inc[u]=[]; g.rev[u]=u
    for a,b in [('A','B'),('B','D'),('A','C'),('C','D')]:
        g.out[a].append(b); g.inc[b].append(a)
        g.edge[(a,b)]=gp.EdgeEvidence(2,0,1)
    return g


def anchor(start: int, end: int, uid: str) -> s20.ExactAnchor:
    return s20.ExactAnchor(start,end,uid)


def test_enumerates_two_bubble_paths() -> None:
    g=bubble_graph()
    paths=s20.enumerate_bounded_paths(g,'A','D',max_edges=3,max_internal_bp=100)
    assert sorted(paths)==sorted([['A','B','D'],['A','C','D']])


def test_k19_internal_evidence_selects_branch() -> None:
    g=bubble_graph()
    paths=s20.enumerate_bounded_paths(g,'A','D',max_edges=3,max_internal_bp=100)
    seq=g.seqs['B']
    choice=s20.choose_path_by_read(
        seq,
        anchor(0,0,'A'),
        anchor(35,35,'D'),
        paths,
        g,
    )
    assert choice is not None
    assert choice.path==('A','B','D')
    assert choice.evidence_k==19
    assert choice.best_hits >= choice.second_hits + 2
    assert choice.source_pos==0 and choice.target_pos==35


def test_flank_only_evidence_does_not_choose_branch() -> None:
    g=bubble_graph()
    paths=s20.enumerate_bounded_paths(g,'A','D',max_edges=3,max_internal_bp=100)
    seq=g.seqs['A'][-35:]+g.seqs['D'][:35]
    choice=s20.choose_path_by_read(
        seq,
        anchor(0,0,'A'),
        anchor(35,35,'D'),
        paths,
        g,
    )
    assert choice is None


def test_missing_physical_edge_blocks_choice() -> None:
    g=bubble_graph()
    g.edge[('B','D')]=gp.EdgeEvidence(0,0,0)
    paths=s20.enumerate_bounded_paths(g,'A','D',max_edges=3,max_internal_bp=100)
    choice=s20.choose_path_by_read(
        g.seqs['B'],
        anchor(0,0,'A'),
        anchor(35,35,'D'),
        paths,
        g,
    )
    assert choice is None


def test_gap_uses_left_end_and_right_start() -> None:
    g=bubble_graph()
    paths=s20.enumerate_bounded_paths(g,'A','D',max_edges=3,max_internal_bp=100)
    seq=g.seqs['B']
    choice=s20.choose_path_by_read(
        seq,
        anchor(0,5,'A'),
        anchor(35,50,'D'),
        paths,
        g,
    )
    assert choice is not None
    assert choice.source_pos==5
    assert choice.target_pos==35


def main() -> None:
    tests=[
        test_enumerates_two_bubble_paths,
        test_k19_internal_evidence_selects_branch,
        test_flank_only_evidence_does_not_choose_branch,
        test_missing_physical_edge_blocks_choice,
        test_gap_uses_left_end_and_right_start,
    ]
    for test in tests:
        test(); print('PASS',test.__name__)


if __name__=='__main__': main()
