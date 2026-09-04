#!/usr/bin/env python3
"""Biological-identity-only wrapper for BridgeBin's same-genome pair head.

Coverage, TNF/composition, GC, and physical linkage are valuable assembly/binning
signals, but they are intentionally excluded from this learned probability head. They
remain available to the Rust graph/core layers where they can support or veto topology
without contaminating the calibration of ``P(same genome | biology)`` across samples.

The feature schema is unchanged for model-file compatibility: pair-level modality values
and presence bits are set to zero before training or scoring. DNA/gene/architecture/
protein/repertoire/taxonomy features are untouched.
"""

from __future__ import annotations

from typing import Dict, List

import bridgebin_pair_head as base


_original_make_features = base.make_features


def biological_features(
    left: base.ContigFeature, right: base.ContigFeature, pair_row: Dict[str, str]
) -> List[float]:
    values = _original_make_features(left, right, pair_row)
    width = len(base.ALL_MODALITIES)
    for name in base.PAIR_FEATURES:
        index = base.ALL_MODALITIES.index(name)
        values[index] = 0.0
        values[width + index] = 0.0
    return values


base.make_features = biological_features


if __name__ == "__main__":
    raise SystemExit(base.main())
