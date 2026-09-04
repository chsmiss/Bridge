#!/usr/bin/env python3
"""DNABERT-S embedding adapter that avoids truncated terminal windows.

The base DNA embedder advances by stride and can emit a short terminal piece. When only a
small fixed number of windows is sampled, that short tail can become one of the selected
identity windows and create length-dependent representation drift. This adapter keeps the
same CLI/runtime but anchors the final window at ``len(sequence) - window_bp`` so every
window from a contig longer than ``window_bp`` has the same length.

It is kept as an explicit ablation wrapper until the behavior is validated end-to-end.
"""

from __future__ import annotations

from typing import Iterator, Tuple

import bridgebin_dna_embed as base


def full_length_windows(
    sequence: str, window_bp: int, stride_bp: int, min_window_bp: int
) -> Iterator[Tuple[int, str]]:
    if len(sequence) <= window_bp:
        if len(sequence) >= min_window_bp:
            yield 0, sequence
        return

    last_start = len(sequence) - window_bp
    starts = list(range(0, last_start + 1, stride_bp))
    if not starts or starts[-1] != last_start:
        starts.append(last_start)
    for serial, start in enumerate(starts):
        yield serial, sequence[start : start + window_bp]


base.windows = full_length_windows


if __name__ == "__main__":
    raise SystemExit(base.main())
