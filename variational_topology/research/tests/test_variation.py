"""Numpy/SciPy checks of the relation functional's first variation and allocations.

Run: python tests/test_variation.py
"""
import sys
import unittest
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from reflection_queries import discordance, foreground_components
from relation_variation import (analog_gain, bpsk_allocation, components, critical_loads, disagreement,
                                grid_neighbors, influence_map, merge_levels, q_function)


def generic_field(rng, side):
    n = side * side
    return (rng.permutation(n).reshape(side, side) + 0.5) / n, 1.0 / n


class CriticalLoads(unittest.TestCase):
    def test_loads_sum_to_one_and_match_enumeration(self):
        rng = np.random.default_rng(0)
        for _ in range(30):
            f, _ = generic_field(rng, 5)
            load = critical_loads(f)
            self.assertAlmostEqual(load.sum(), 1.0, places=12)
            self.assertGreaterEqual(load.min(), 1.0 / f.size**2 - 1e-15)       # the pair (v, v)
            m = merge_levels(f)
            order = np.argsort(f.ravel())
            owner = order[np.searchsorted(f.ravel()[order], m)]
            brute = np.bincount(owner.ravel(), minlength=f.size) / f.size**2
            np.testing.assert_allclose(brute, load.ravel(), atol=1e-14)

    def test_local_law_is_a_weighted_l1_norm(self):
        rng = np.random.default_rng(1)
        for _ in range(30):
            f, gap = generic_field(rng, 5)
            load, m = critical_loads(f), merge_levels(f)
            delta = rng.uniform(-0.499 * gap, 0.499 * gap, f.shape)
            direct = np.abs(merge_levels(f + delta) - m).mean()
            self.assertAlmostEqual(direct, (load * np.abs(delta)).sum(), places=14)
            cdf = lambda t: t**2                                               # another threshold law on [0, 1]
            direct = np.abs(cdf(merge_levels(f + delta)) - cdf(m)).mean()
            self.assertAlmostEqual(direct, (load * np.abs(cdf(f + delta) - cdf(f))).sum(), places=14)

    def test_law_can_fail_beyond_half_the_gap(self):
        rng = np.random.default_rng(2)
        failures = 0
        for _ in range(30):
            f, gap = generic_field(rng, 5)
            delta = rng.uniform(-3 * gap, 3 * gap, f.shape)
            direct = np.abs(merge_levels(f + delta) - merge_levels(f)).mean()
            failures += abs(direct - (critical_loads(f) * np.abs(delta)).sum()) > 1e-9
        self.assertGreater(failures, 0)


def simple_paths(nbrs, a, b):
    stack = [(a, (a,))]
    while stack:
        node, path = stack.pop()
        if node == b:
            yield path
            continue
        for nxt in nbrs[node]:
            if nxt not in path:
                stack.append((nxt, path + (nxt,)))


class DirectionalDerivative(unittest.TestCase):
    def test_minimax_derivative_with_ties(self):
        rng = np.random.default_rng(3)
        nbrs = grid_neighbors(3, 3)
        for _ in range(8):
            f = rng.choice([0.0, 0.5, 1.0], size=9)
            delta = rng.normal(size=9)
            step = 1e-6
            m0 = merge_levels(f.reshape(3, 3))
            m1 = merge_levels((f + step * delta).reshape(3, 3))
            for a, b in product(range(9), repeat=2):
                paths = [(a,)] if a == b else list(simple_paths(nbrs, a, b))
                best = min(max(f[list(p)]) for p in paths)
                formula = min(max(delta[v] for v in p if f[v] == best)
                              for p in paths if max(f[list(p)]) == best)
                self.assertAlmostEqual((m1[a, b] - m0[a, b]) / step, formula, places=6)


class BinaryInfluence(unittest.TestCase):
    def test_fast_count_matches_reference(self):
        rng = np.random.default_rng(4)
        for _ in range(50):
            x, y = rng.random((6, 6)) < 0.5, rng.random((6, 6)) < 0.5
            ref = discordance(foreground_components(x), foreground_components(y), np.ones(x.shape, bool))
            self.assertAlmostEqual(disagreement(components(x), components(y)), ref, places=15)

    def test_first_order_remainder_bounds(self):
        rng = np.random.default_rng(5)
        for _ in range(3):
            x = rng.random((3, 4)) < 0.5
            base, inf = components(x), influence_map(x)
            eps = rng.uniform(0.0, 0.03, x.shape)
            total = 0.0
            for bits in product((False, True), repeat=x.size):
                z = np.array(bits).reshape(x.shape)
                prob = np.prod(np.where(z, eps, 1 - eps))
                total += prob * disagreement(base, components(x ^ z))
            first = (eps * inf).sum()
            s = eps.sum()
            self.assertLessEqual(total - first, s * s / 2 + 1e-15)
            self.assertGreaterEqual(total - first, -s * first - 1e-15)
            self.assertGreaterEqual(inf.min(), 1.0 / x.size**2 - 1e-15)


class Allocation(unittest.TestCase):
    def test_analog_rule_and_gain(self):
        rng = np.random.default_rng(6)
        for _ in range(10):
            f, _ = generic_field(rng, 5)
            c = critical_loads(f).ravel()
            n, power = c.size, 7.0
            best = power * c ** (2 / 3) / (c ** (2 / 3)).sum()
            value = (c / np.sqrt(best)).sum()
            uniform = (c / np.sqrt(power / n)).sum()
            self.assertAlmostEqual(value / uniform, analog_gain(c), places=12)
            self.assertLessEqual(analog_gain(c), 1.0)
            for _ in range(200):
                p = power * rng.dirichlet(np.ones(n))
                self.assertGreaterEqual((c / np.sqrt(p)).sum(), value - 1e-9)
        self.assertAlmostEqual(analog_gain(np.full(16, 1 / 16)), 1.0, places=12)

    def test_bpsk_allocation_is_stationary_and_optimal(self):
        rng = np.random.default_rng(7)
        for snr_db in (2.0, 6.0):
            w = rng.uniform(1e-6, 3e-3, 40) * rng.choice([1.0, 1e-3], 40)
            snr = 10 ** (snr_db / 10)
            p = bpsk_allocation(w, snr)
            self.assertAlmostEqual(p.mean(), snr, places=9)
            slope = w * np.exp(-p) / (2 * np.sqrt(np.pi * p))                 # -d/dp of w Q(sqrt(2p))
            np.testing.assert_allclose(slope, slope[0], rtol=1e-8)
            value = (w * q_function(np.sqrt(2 * p))).sum()
            self.assertLessEqual(value, (w * q_function(np.sqrt(2 * snr))).sum())
            for _ in range(300):
                trial = np.clip(p * np.exp(0.2 * rng.normal(size=p.size)), 1e-9, None)
                trial *= snr / trial.mean()
                self.assertGreaterEqual((w * q_function(np.sqrt(2 * trial))).sum(), value - 1e-15)


if __name__ == "__main__":
    unittest.main()
