"""Regression checks for mathematical edge cases and independent libraries."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code/topology'), str(ROOT / 'code/embeddings')]
import tda
import sync
import reeb
import run_topology
import topology_exp


class TopologyTests(unittest.TestCase):
    def test_bottleneck_against_independent_library(self):
        import persim
        rng = np.random.default_rng(17)
        for _ in range(40):
            pair = []
            for _ in range(2):
                n = int(rng.integers(0, 10))
                birth = rng.uniform(0, 2, n)
                pair.append(np.c_[birth, birth + rng.exponential(.4, n)])
            self.assertAlmostEqual(tda.bottleneck(*pair), persim.bottleneck(*pair), places=12)

    def test_essential_classes(self):
        self.assertEqual(tda.bottleneck([[0, np.inf]], []), np.inf)
        self.assertEqual(tda.bottleneck([[0, np.inf]], [[.2, np.inf]]), .2)

    def test_duplicate_noise_reference(self):
        D = np.zeros((8, 8))
        expected = np.ones(7) / np.sqrt(2)
        ref = tda.noise_reference(D, 1)
        np.testing.assert_allclose(tda.h0_deaths(D=ref), expected)
        transformed = tda.phi_reparam(tda.h0_diagram(tda.h0_deaths(D=D)), 1, 0)
        np.testing.assert_allclose(transformed[:-1, 1], expected)

    def test_zero_mst_edges_are_retained(self):
        np.testing.assert_allclose(tda.h0_deaths(X=np.array([[0.], [0.], [2.]])), [1., 0.])
        self.assertEqual(len(tda.h0_deaths(X=np.zeros((5, 2)))), 4)

    def test_native_alpha_against_gudhi(self):
        rng = np.random.default_rng(18)
        for n in [4, 13, 50]:
            for _ in range(3):
                X = rng.normal(size=(n, 2))
                actual = tda.alpha2d(X, backend='native')
                expected = tda.alpha2d(X, backend='gudhi')
                self.assertLess(max(tda.bottleneck_by_degree(actual, expected)), 2e-10)

    def test_collapsed_interval_has_no_loop(self):
        X = np.c_[np.linspace(-1, 1, 30), np.zeros(30)]
        d = tda.alpha2d(X)
        self.assertEqual(len(d[1]), 0)
        self.assertEqual(tda.n_significant(d[0], .1), 1)

    def test_maxmin_unique_on_duplicates(self):
        ids, radius = tda.maxmin(np.zeros((8, 2)), 8)
        self.assertEqual(len(np.unique(ids)), 8)
        self.assertEqual(radius, 0)

    def test_mixed_radix_nonpower_of_two(self):
        reference = np.array([[1., 2.], [3., 4.]])
        for offsets in [np.full((2, 2), -.1), np.full((2, 2), .1), np.array([[-.07, .01], [.03, .09]])]:
            msg, m = sync.encode(reference, reference + offsets, .1, .04)
            self.assertEqual(m, 3)
            self.assertEqual(len(msg), 7)  # ceil(log2(3^4)), not 4 ceil(log2(3)).
            recovered = sync.decode(reference, msg, .1, m)
            self.assertLessEqual(tda.bottleneck(reference + offsets, recovered), .04)

    def test_unused_codewords_rejected(self):
        with self.assertRaises(ValueError):
            sync.decode([[1, 2]], '1111', .1, 3)

    def test_identity_alignment_reorders_caches(self):
        raw = {'A': {'ids': np.array(['x','y','z','w']), 'X': np.array([[0],[1],[2],[3]]), 'labels': np.array([0,1,0,1])},
               'B': {'ids': np.array(['z','x','w','y']), 'X': np.array([[20],[0],[30],[10]]), 'labels': np.array([0,0,1,1])}}
        ids, X, y = run_topology.align(raw)
        np.testing.assert_array_equal(X['B'], X['A'] * 10)

    def test_label_disagreement_rejected(self):
        raw = {'A': {'ids': np.array(['a','b','c','d']), 'X': np.zeros((4,2)), 'labels': np.array([0,0,1,1])},
               'B': {'ids': np.array(['a','b','c','d']), 'X': np.zeros((4,2)), 'labels': np.array([1,0,1,1])}}
        with self.assertRaises(ValueError):
            run_topology.align(raw)

    def test_missing_ids_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'm__d'; path.mkdir()
            np.savez(path / 'arrays.npz', image=np.zeros((4,2)))
            with self.assertRaisesRegex(ValueError, 'missing item_ids'):
                run_topology.load(Path(tmp), 'm', 'd', 'image')

    def test_probe_requires_reference(self):
        with self.assertRaises(ValueError):
            topology_exp.probe_score(np.eye(3), np.array([0,1,0]))

    def test_mapper_requires_graph_cover(self):
        with self.assertRaises(ValueError):
            reeb.mapper_graph(np.eye(3), np.arange(3), 3, .8, .5)


if __name__ == '__main__':
    unittest.main(verbosity=2)
