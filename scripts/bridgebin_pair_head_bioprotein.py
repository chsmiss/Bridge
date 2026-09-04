#!/usr/bin/env python3
"""Biological pair head with ESM-C set-to-set protein evidence.

Coverage/composition/GC/physical linkage remain outside the learned biological identity
head.  When candidate rows contain ``protein_set_similarity`` this wrapper injects that
set-level score into the existing ``protein`` modality, overriding the less informative
pooled-contig protein cosine.  The model schema stays unchanged.
"""

from __future__ import annotations

from typing import Dict, List

import bridgebin_pair_head as base


_original_make_features = base.make_features


def biological_protein_features(
    left: base.ContigFeature, right: base.ContigFeature, pair_row: Dict[str, str]
) -> List[float]:
    values = _original_make_features(left, right, pair_row)
    width = len(base.ALL_MODALITIES)
    for name in base.PAIR_FEATURES:
        index = base.ALL_MODALITIES.index(name)
        values[index] = 0.0
        values[width + index] = 0.0

    protein_score = base.optional_pair_value(
        pair_row, ("protein_set_similarity", "protein_pair_similarity")
    )
    if protein_score is not None:
        confidence = base.optional_pair_value(
            pair_row, ("protein_confidence", "protein_set_confidence")
        )
        if confidence is None:
            confidence = min(left.protein_confidence, right.protein_confidence)
        confidence = max(0.0, min(1.0, confidence))
        index = base.ALL_MODALITIES.index("protein")
        values[index] = protein_score * confidence
        values[width + index] = confidence
    return values


base.make_features = biological_protein_features


if __name__ == "__main__":
    raise SystemExit(base.main())
