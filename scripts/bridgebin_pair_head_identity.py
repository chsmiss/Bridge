#!/usr/bin/env python3
"""DNABERT-S identity-head wrapper with geometry and uncertainty kept separate.

The generic BridgeBin pair head multiplies a modality similarity by its confidence before
standardization.  That is sensible for missing/noisy generic modalities, but for genome
identity it can reorder pairs merely because one contig has a shorter terminal window.
This wrapper keeps DNA cosine as the learned identity coordinate and turns DNA confidence
into a binary presence bit for the probability head.  Coverage, composition, GC and
physical-link features remain excluded exactly as in the Biological Brain firewall.

The input/output model schema is unchanged, so the normal training/scoring implementation
and threshold metadata remain reusable. Models trained with this wrapper must also be
scored with this wrapper.
"""

from __future__ import annotations

from typing import Dict, List

import bridgebin_pair_head as base


_original_make_features = base.make_features


def identity_features(
    left: base.ContigFeature, right: base.ContigFeature, pair_row: Dict[str, str]
) -> List[float]:
    values = _original_make_features(left, right, pair_row)
    width = len(base.ALL_MODALITIES)

    # Hard modality firewall: cheap graph evidence never enters the learned biological
    # identity probability.
    for name in base.PAIR_FEATURES:
        index = base.ALL_MODALITIES.index(name)
        values[index] = 0.0
        values[width + index] = 0.0

    # The generic builder stores DNA as cosine * confidence and confidence separately.
    # Recover the raw cosine and use the second slot as presence only. Reliability can be
    # consumed by downstream gating without changing the identity ranking itself.
    dna_index = base.ALL_MODALITIES.index("dna")
    dna_confidence = values[width + dna_index]
    if dna_confidence > 0.0:
        values[dna_index] = max(0.0, min(1.0, values[dna_index] / dna_confidence))
        values[width + dna_index] = 1.0
    else:
        values[dna_index] = 0.0
        values[width + dna_index] = 0.0

    return values


base.make_features = identity_features


if __name__ == "__main__":
    raise SystemExit(base.main())
