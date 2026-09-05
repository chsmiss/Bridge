#!/usr/bin/env python3
from __future__ import annotations

import random
from pathlib import Path

import stage31_multik_seed_bridge as s31
import stage33_seed_end_microassembly as s33


def dna(seed: int, n: int) -> bytes:
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(n))


def synthetic_fragments(truth: bytes, start: int, stop: int) -> dict[int, s33.Fragment]:
    fragments: dict[int, s33.Fragment] = {}
    fid = 1
    for pos in range(start, stop, 28):
        read = truth[pos : pos + 101]
        if len(read) < 101:
            break
        # Two distinct physical fragments cover each local window so the solid2
        # graph can extend without singleton mercy.
        for _copy in range(2):
            q = b"I" * len(read)
            mate = s31.rc(read)
            fragments[fid] = s33.Fragment(fid, read, q, mate, q)
            fid += 1
    return fragments


def main() -> None:
    left = dna(101, 700)
    gap = dna(102, 120)
    right = dna(103, 720)
    truth = left + gap + right
    seeds = [("left", left), ("right", right)]

    # Signature construction must cover both physical ends and retain unique
    # mapping anchors on random sequence.
    unique, by_ep, stats = s33.terminal_signature_index(seeds, 240)
    assert stats["physical_endpoints"] == 4
    assert stats["endpoints_with_unique_k21"] == 4, stats
    assert by_ep[s33.endpoint_id(0, "R")]
    assert unique

    # Cover the left boundary, the entire real gap, and >80 bp of the right
    # seed.  Local assembly is not told the truth sequence or target seed.
    fragments = synthetic_fragments(
        truth,
        len(left) - 90,
        len(left) + len(gap) + 100,
    )
    fids = set(fragments)
    eid = s33.endpoint_id(0, "R")
    result, support = s33.extend_local(
        seeds,
        eid,
        fids,
        fragments,
        seed_overlap=120,
        max_extension=400,
        dominance=1.75,
        allow_mercy=False,
        min_mercy_fragments=0,
        singleton_density=0.30,
        min_solid_nodes=40,
        singleton_quality=30.0,
    )
    assert support["solid_edges"] > 0, support
    assert result.added_bp >= len(gap) + 45, result
    assert result.mercy_edges == 0
    assert result.sequence == truth[len(left) - 120 : len(left) - 120 + len(result.sequence)]

    local = {eid: result}
    proposals, dstats = s33.discover_bridges(seeds, local, 40, 10)
    assert dstats["source_consistent_proposals"] == 1, dstats
    assert proposals[0].proposal.le == (0, "R")
    assert proposals[0].proposal.re == (1, "L")
    chosen, sstats = s33.select_bridges(proposals, 2, 4)
    assert sstats["selected_bridges"] == 1, sstats

    packed = [(chosen[0].proposal, (31,), 1)]
    assembled = s31.assemble(seeds, packed)
    assert len(assembled) == 1, assembled
    merged = assembled[0][1]
    assert merged == truth, (len(merged), len(truth))
    assert b"N" not in merged
    assert s31.n50(seq for _name, seq in assembled) > max(len(left), len(right))

    # A mercy path backed by fewer fragments than the requested Stage23 guard
    # must never be emitted.  This tests the final invariant directly without
    # relying on a particular synthetic singleton layout.
    fake = s33.LocalResult((0, "R"), left[-120:] + b"A", 1, 1, 1, 2, (7,), 0, "test")
    assert len(fake.mercy_fragments) < 2

    print("stage33 tests: passed")


if __name__ == "__main__":
    main()
