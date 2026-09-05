#!/usr/bin/env python3
from __future__ import annotations

import random
import tempfile
from collections import Counter
from pathlib import Path

import graph_path_phaser as gp
import stage31_multik_seed_bridge as s31
import stage32_targeted_gfa_seed_bridge as s32


def dna(seed: int, n: int) -> bytes:
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(n))


def main() -> None:
    left = dna(11, 700)
    right = dna(12, 750)
    middle = dna(13, 420)
    core = left[-180:] + middle + right[:180]
    rev = s31.rc(core)
    with tempfile.TemporaryDirectory() as td:
        gfa = Path(td) / "graph.gfa"
        gfa.write_text(
            "H\tVN:Z:1.0\n"
            + f"S\tf\t{core.decode()}\tKC:f:10\n"
            + f"S\tr\t{rev.decode()}\tKC:f:10\n"
        )
        graph = gp.Graph.from_gfa(gfa)
        index = gp.KmerIndex(graph, 31)
        seeds = [("left", left), ("right", right)]
        raw = s32.anchor_candidates(seeds, graph, index, 80, 180)
        anchors, astats = s32.choose_anchors(raw)
        assert astats["chosen_anchors"] >= 4, astats

        all_props = []
        for name in s32.PASS_ORDER:
            vals, stats = s32.discover_pass(
                seeds,
                graph,
                anchors,
                Counter(),
                Counter(),
                Counter(),
                name,
                8,
                2000,
            )
            assert stats["proposals"] >= 1, (name, stats)
            all_props.extend(vals)

        candidates, agg = s32.aggregate_edges(all_props, 2, True)
        assert agg["candidate_edges"] >= 1, agg
        chosen, sel = s32.select_edges(candidates, 2, 4)
        assert sel["selected_bridges"] == 1, sel
        packed = [
            (edge.p, tuple(31 for _ in edge.passes), edge.evidence_rank[0])
            for edge in chosen
        ]
        assembled = s31.assemble(seeds, packed)
        assert len(assembled) == 1, assembled
        seq = assembled[0][1]
        assert left in seq and right in seq and middle in seq
        assert b"N" not in seq
        assert s31.n50(x for _, x in assembled) > max(len(left), len(right))

    print("stage32 tests: passed")


if __name__ == "__main__":
    main()
