#!/usr/bin/env python3
from __future__ import annotations

from seed_locked_extensions import seed_locked_extensions


def dna(length: int, offset: int = 0) -> bytes:
    # Deterministic non-periodic-enough sequence for exact-overlap unit tests.
    x = 0x12345678 + offset * 0x9E3779B1
    out = bytearray()
    alphabet = b"ACGT"
    for _ in range(length):
        x ^= (x << 13) & 0xFFFFFFFF
        x ^= x >> 17
        x ^= (x << 5) & 0xFFFFFFFF
        out.append(alphabet[x & 3])
    return bytes(out)


def test_unique_right_extension() -> None:
    seed = dna(320, 1)
    extra = dna(90, 2)
    candidate = seed[-160:] + extra
    output, stats = seed_locked_extensions(
        [("seed", seed)],
        [("candidate", candidate)],
        min_overlap=120,
        overlap_margin=20,
        seed_length=31,
        min_extension=20,
        max_seed_occurrences=64,
    )
    assert len(output) == 1
    assert output[0][1] == seed + extra
    assert stats["accepted_extensions"] == 1
    assert stats["added_bp"] == len(extra)
    assert stats["output_n50"] == len(seed) + len(extra)


def test_seed_is_not_deleted_without_extension() -> None:
    seed = dna(320, 3)
    contained = seed[40:280]
    output, stats = seed_locked_extensions(
        [("seed", seed)],
        [("candidate", contained)],
        min_overlap=120,
        overlap_margin=20,
        seed_length=31,
        min_extension=20,
        max_seed_occurrences=64,
    )
    assert output[0][1] == seed
    assert stats["accepted_extensions"] == 0
    assert stats["added_bp"] == 0


def test_two_seed_bridge_is_rejected() -> None:
    left = dna(320, 4)
    right = dna(330, 5)
    bridge = left[-150:] + dna(80, 6) + right[:150]
    output, stats = seed_locked_extensions(
        [("left", left), ("right", right)],
        [("bridge", bridge)],
        min_overlap=120,
        overlap_margin=20,
        seed_length=31,
        min_extension=20,
        max_seed_occurrences=64,
    )
    assert [sequence for _header, sequence in output] == [left, right]
    assert stats["ambiguous_bridge_candidates_rejected"] == 1
    assert stats["accepted_extensions"] == 0


def test_both_seed_ends_can_extend_with_distinct_candidates() -> None:
    seed = dna(360, 7)
    left_extra = dna(70, 8)
    right_extra = dna(75, 9)
    left_candidate = left_extra + seed[:150]
    right_candidate = seed[-150:] + right_extra
    output, stats = seed_locked_extensions(
        [("seed", seed)],
        [("left", left_candidate), ("right", right_candidate)],
        min_overlap=120,
        overlap_margin=20,
        seed_length=31,
        min_extension=20,
        max_seed_occurrences=64,
    )
    assert output[0][1] == left_extra + seed + right_extra
    assert stats["accepted_extensions"] == 2
    assert stats["left_extensions"] == 1
    assert stats["right_extensions"] == 1


if __name__ == "__main__":
    test_unique_right_extension()
    test_seed_is_not_deleted_without_extension()
    test_two_seed_bridge_is_rejected()
    test_both_seed_ends_can_extend_with_distinct_candidates()
    print("seed-locked extension tests: passed")
