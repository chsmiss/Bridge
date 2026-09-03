#!/usr/bin/env python3
from __future__ import annotations

import unittest
from collections import Counter

import stage15_structural_recovery as s15


class Stage15StructuralRecoveryTests(unittest.TestCase):
    def test_soft_context_increment_uses_excess_not_double_count(self) -> None:
        baseline = Counter({("a", "b"): 3, ("x", "y", "z"): 2})
        k21 = Counter({("a", "b"): 5, ("x", "y", "z"): 3, ("c", "d"): 1})
        k17 = Counter({("a", "b"): 4, ("x", "y", "z"): 6, ("e", "f"): 2})
        result = s15.soft_context_increment(baseline, [k21, k17])
        self.assertEqual(result[("a", "b")], 2)
        self.assertEqual(result[("x", "y", "z")], 4)
        self.assertEqual(result[("e", "f")], 2)
        self.assertNotIn(("c", "d"), result)

    def test_soft_context_increment_caps_weight(self) -> None:
        result = s15.soft_context_increment(
            Counter(), [Counter({("a", "b", "c"): 99})], max_weight=3
        )
        self.assertEqual(result[("a", "b", "c")], 3)

    def test_soft_context_increment_ignores_singletons(self) -> None:
        result = s15.soft_context_increment(Counter(), [Counter({("a",): 10})])
        self.assertFalse(result)


if __name__ == "__main__":
    unittest.main()
