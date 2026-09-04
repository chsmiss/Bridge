#!/usr/bin/env python3
"""OOD-safe entry point for BridgeBin's multimodal same-genome pair head.

This keeps the v3 pair-head feature schema and calibration logic, but replaces the
unbounded z-score transform with a conservative normalization suitable for transferring
between communities/libraries:

* similarity features use a minimum scale floor so near-constant training modalities
  cannot explode at inference time;
* presence/confidence features use a larger floor;
* all standardized coordinates are clipped before entering the linear head.

The wrapper intentionally reuses ``bridgebin_pair_head`` for parsing, training,
calibration, and scoring so model files remain compatible with the production pipeline.
"""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple

import bridgebin_pair_head as base


SIMILARITY_SCALE_FLOOR = 0.05
PRESENCE_SCALE_FLOOR = 0.10
Z_CLIP = 4.0


def robust_fit_standardizer(examples: Sequence[base.Example]) -> Tuple[List[float], List[float]]:
    width = len(base.FEATURE_NAMES)
    means = [0.0] * width
    for example in examples:
        for index, value in enumerate(example.features):
            means[index] += value
    means = [value / len(examples) for value in means]

    variances = [0.0] * width
    for example in examples:
        for index, value in enumerate(example.features):
            variances[index] += (value - means[index]) ** 2

    scales: List[float] = []
    similarity_width = len(base.ALL_MODALITIES)
    for index, variance in enumerate(variances):
        observed = math.sqrt(variance / len(examples))
        floor = SIMILARITY_SCALE_FLOOR if index < similarity_width else PRESENCE_SCALE_FLOOR
        scales.append(max(observed, floor))
    return means, scales


def robust_standardized(
    values: Sequence[float], means: Sequence[float], scales: Sequence[float]
) -> List[float]:
    result: List[float] = []
    for value, mean, scale in zip(values, means, scales):
        z = (value - mean) / max(scale, 1e-12)
        result.append(max(-Z_CLIP, min(Z_CLIP, z)))
    return result


base.fit_standardizer = robust_fit_standardizer
base.standardized = robust_standardized


if __name__ == "__main__":
    raise SystemExit(base.main())
