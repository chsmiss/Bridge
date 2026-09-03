#!/usr/bin/env python3
from __future__ import annotations

import unittest
from collections import Counter, defaultdict
from pathlib import Path

import graph_path_phaser as gp
import stage14_amplified_methods as s14


def graph_from_edges(seqs: dict[str, str], edges: list[tuple[str, str]], coverage=None, k=3):
    g = gp.Graph()
    g.k = k
    g.seqs = dict(seqs)
    g.coverage = {uid: float((coverage or {}).get(uid, 10.0)) for uid in seqs}
    g.out = defaultdict(list)
    g.inc = defaultdict(list)
    g.edge = {}
    g.rev = {uid: uid for uid in seqs}
    for a, b in edges:
        g.out[a].append(b)
        g.inc[b].append(a)
        g.edge[(a, b)] = gp.EdgeEvidence(direct=4, gapped=0, pairs=1)
    for uid in seqs:
        g.out[uid] = sorted(set(g.out.get(uid, [])))
        g.inc[uid] = sorted(set(g.inc.get(uid, [])))
    return g


class Stage14Tests(unittest.TestCase):
    def test_unique_bounded_bridge_accepts_only_unique_path(self):
        g = graph_from_edges(
            {x: 'A' * 20 for x in ('s', 'a', 'b', 't')},
            [('s', 'a'), ('a', 't')],
        )
        self.assertEqual(s14.unique_bounded_bridge(g, 's', 't'), ['s', 'a', 't'])
        g.out['s'].append('b')
        g.inc['t'].append('b')
        g.edge[('s', 'b')] = gp.EdgeEvidence(direct=4, pairs=1)
        g.edge[('b', 't')] = gp.EdgeEvidence(direct=4, pairs=1)
        g.out['b'] = ['t']
        self.assertIsNone(s14.unique_bounded_bridge(g, 's', 't'))

    def test_extend_threaded_path_walks_unique_chain(self):
        g = graph_from_edges(
            {x: ('ACGT' * 8)[:28] for x in ('a', 'b', 'c', 'd')},
            [('a', 'b'), ('b', 'c'), ('c', 'd')],
        )
        extended = s14.extend_threaded_path(g, ['b', 'c'], Counter(), max_nodes=10, max_bp=500)
        self.assertEqual(extended, ['a', 'b', 'c', 'd'])

    def test_context_extension_refuses_ambiguous_unsupported_branch(self):
        g = graph_from_edges(
            {x: ('ACGT' * 8)[:28] for x in ('s', 'a', 'b')},
            [('s', 'a'), ('s', 'b')],
        )
        g.edge[('s', 'a')] = gp.EdgeEvidence()
        g.edge[('s', 'b')] = gp.EdgeEvidence()
        self.assertIsNone(s14._choose_context_extension(g, ['s'], Counter(), forward=True))
        ctx = Counter({('s', 'a'): 4, ('s', 'b'): 1})
        self.assertIsNone(s14._choose_context_extension(g, ['s'], ctx, forward=True))

    def test_variable_k_promotes_branch_reduction(self):
        base = s14.LocalAssembly(1, 17, Path('/x'), 1000, 200, 10, 20)
        high = s14.LocalAssembly(1, 25, Path('/y'), 750, 260, 6, 5)
        chosen, reason = s14.choose_variable_k([base, high])
        self.assertIsNotNone(chosen)
        assert chosen is not None
        self.assertEqual(chosen.k, 25)
        self.assertEqual(reason, 'promote_higher_k')

    def test_variable_k_keeps_base_when_high_k_loses_too_much(self):
        base = s14.LocalAssembly(1, 17, Path('/x'), 1000, 200, 10, 20)
        high = s14.LocalAssembly(1, 31, Path('/y'), 300, 500, 2, 1)
        chosen, reason = s14.choose_variable_k([base, high])
        self.assertIsNotNone(chosen)
        assert chosen is not None
        self.assertEqual(chosen.k, 17)
        self.assertEqual(reason, 'keep_k17')

    def test_long_tangle_prefers_far_reconvergence(self):
        seqs = {x: ('ACGT' * 30)[:100] for x in ('s','a','b','m','c','d','e','f','z')}
        edges = [
            ('s','a'),('s','b'),('a','m'),('b','m'),
            ('m','c'),('m','d'),('c','e'),('d','f'),('e','z'),('f','z')
        ]
        g = graph_from_edges(seqs, edges, k=31)
        tangles = s14.discover_long_tangles(
            g, min_bp=300, max_bp=1000, max_depth=12, max_paths=8, max_tangles=4
        )
        self.assertTrue(tangles)
        target = next((t for t in tangles if t.source == 's' and t.sink == 'z'), None)
        self.assertIsNotNone(target)
        assert target is not None
        self.assertGreaterEqual(target.min_bp, 300)

    def test_thread_fraction_measures_path_coverage(self):
        ctx = Counter({('a','b'): 3, ('b','c'): 2})
        supported, frac, maximum = s14.path_thread_fraction(ctx, ['a','b','c','d'])
        self.assertEqual(supported, 2)
        self.assertAlmostEqual(frac, 2/3)
        self.assertEqual(maximum, 3)

    def test_seed_assignment_requires_clear_margin(self):
        self.assertEqual(s14.assign_seed(Counter({1: 4, 2: 1}), min_hits=2, margin=2), 1)
        self.assertIsNone(s14.assign_seed(Counter({1: 3, 2: 2}), min_hits=2, margin=2))


if __name__ == '__main__':
    unittest.main()

# Trigger Stage14 benchmark through the contents API.
