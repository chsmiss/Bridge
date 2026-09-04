#!/usr/bin/env python3
"""BridgeBin Biological Brain v2.2 experimental production wrapper.

Compared with the existing v2.1 orchestrator this wrapper makes two conservative changes
without duplicating the pipeline implementation:

1. expensive DNA-LM inference is restricted to contigs that actually appear in the sparse
   candidate-pair table and uses an even cap on windows per contig;
2. the learned pair probability is biological identity only. Coverage/composition/GC and
   physical linkage stay in the graph/core layers instead of entering the calibrated
   same-genome head.

The wrapper monkey-patches only the orchestration hooks; all normal input validation,
feature extraction, manifest writing, and v2.1 Rust refinement remain in the base module.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import bridgebin_biobrain_pipeline as base


_candidate_pairs: Optional[Path] = None
_original_mine_candidates = base.mine_candidates
_original_run = base.run
_original_score_pairs = base.score_pairs


def mine_candidates(args, assignments: Path, work: Path) -> Path:
    global _candidate_pairs
    _candidate_pairs = _original_mine_candidates(args, assignments, work)
    return _candidate_pairs


def run(command: Sequence[object], log: Optional[Path] = None) -> None:
    argv = list(command)
    if len(argv) >= 2 and str(argv[1]).endswith("bridgebin_dna_embed.py"):
        if _candidate_pairs is not None and "--pairs" not in [str(value) for value in argv]:
            argv.extend(["--pairs", _candidate_pairs])
        if "--max-windows-per-contig" not in [str(value) for value in argv]:
            argv.extend(["--max-windows-per-contig", 2])
    _original_run(argv, log)


def score_pairs(args, features: Path, candidates: Path, work: Path) -> Path:
    output = work / "pair_scores.tsv"
    run(
        [
            base.sys.executable,
            base.SCRIPTS / "bridgebin_pair_head_bio.py",
            "score",
            "--features",
            features,
            "--pairs",
            candidates,
            "--model",
            args.pair_model,
            "--output",
            output,
            "--model-name",
            "bridgebin-biobrain-v22-bio-only",
        ],
        work / "09_pair_head.log",
    )
    return output


base.mine_candidates = mine_candidates
base.run = run
base.score_pairs = score_pairs


if __name__ == "__main__":
    raise SystemExit(base.main())
