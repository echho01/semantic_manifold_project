"""NumPy-only checks of the orbit, merge-function and fading statements (no torch or GUDHI)."""

import itertools
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.special import exp1

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from orbit_rdf import DIHEDRAL, blahut_point, distortion_matrix, orbit, orbit_rate


def merge_matrix(f):
    """Minimax path values m_f(a,b) on the four-neighbor grid by Floyd-Warshall."""
    h, w = f.shape
    n = h * w
    val = f.ravel()
    m = np.full((n, n), np.inf)
    for v in range(n):
        m[v, v] = val[v]
        y, x = divmod(v, w)
        for a, b in ((y + 1, x), (y, x + 1)):
            if a < h and b < w:
                u = a * w + b
                m[v, u] = m[u, v] = max(val[v], val[u])
    for k in range(n):
        m = np.minimum(m, np.maximum(m[:, [k]], m[[k], :]))
    return m


def h0_diagram(f):
    """Sublevel H0 pairs by the elder rule with a union-find over vertices."""
    h, w = f.shape
    val = f.ravel()
    parent, birth, seen, pairs = {}, {}, set(), []

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for v in np.argsort(val, kind="stable"):
        v = int(v)
        parent[v], birth[v] = v, val[v]
        seen.add(v)
        y, x = divmod(v, w)
        for a, b in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
            u = a * w + b
            if 0 <= a < h and 0 <= b < w and u in seen:
                ru, rv = find(u), find(v)
                if ru != rv:
                    young, old = (ru, rv) if birth[ru] > birth[rv] else (rv, ru)
                    if val[v] > birth[young]:
                        pairs.append((birth[young], val[v]))
                    parent[young] = old
    pairs.append((val.min(), np.inf))
    return sorted(pairs)


class Orbit(unittest.TestCase):
    def test_merge_function_and_automorphisms(self):
        rng = np.random.default_rng(7)
        f = rng.uniform(size=(4, 4))
        m = merge_matrix(f)
        n = f.size
        # Strong triangle inequality and the diagonal convention.
        self.assertTrue(np.all(m[:, None, :] <= np.maximum(m[:, :, None], m[None, :, :]) + 1e-12))
        np.testing.assert_allclose(np.diag(m), f.ravel())
        index = np.arange(n).reshape(f.shape)
        for name, sigma in DIHEDRAL.items():
            g = sigma(f)
            perm = sigma(index).ravel()  # g.ravel()[v] = f.ravel()[perm[v]]
            np.testing.assert_allclose(merge_matrix(g), m[np.ix_(perm, perm)], err_msg=name)
            # The unlabelled H0 diagram does not change on the orbit.
            np.testing.assert_allclose(np.array(h0_diagram(g)), np.array(h0_diagram(f)))

    def test_grey_orbit_distortion_is_balanced(self):
        rng = np.random.default_rng(11)
        f = rng.uniform(size=(4, 4))
        members = [np.ascontiguousarray(s(f)) for s in DIHEDRAL.values()]
        merges = [merge_matrix(g) for g in members]
        rho = np.array([[np.abs(a - b).mean() for b in merges] for a in merges])
        rows = np.sort(rho, axis=1)
        np.testing.assert_allclose(rows, np.tile(rows[0], (8, 1)), atol=1e-12)
        for beta in (3.0, 20.0, 80.0):
            d, r = blahut_point(rho, beta)
            self.assertAlmostEqual(r, orbit_rate(rho[0], d), places=7)

    def test_binary_orbit_formula_matches_blahut_and_enumeration(self):
        rng = np.random.default_rng(3)
        for _ in range(20):
            mask = rng.random((5, 5)) < 0.45
            members = orbit(mask)
            rho = distortion_matrix(members)
            # Brute-force enumeration through the merge function of the 0/1 score.
            merges = [merge_matrix(1.0 - m) for m in members]
            brute = np.array([[np.mean((a < 0.5) != (b < 0.5)) for b in merges] for a in merges])
            np.testing.assert_allclose(rho, brute, atol=1e-12)
            rows = np.sort(rho, axis=1)
            np.testing.assert_allclose(rows, np.tile(rows[0], (len(members), 1)), atol=1e-12)
            for beta in (4.0, 30.0, 200.0):
                d, r = blahut_point(rho, beta)
                self.assertAlmostEqual(r, orbit_rate(rho[0], d), places=7)
            self.assertAlmostEqual(orbit_rate(rho[0], rho[0].mean()), 0.0)

    def test_two_template_special_case(self):
        alpha = 0.3
        for delta in (0.01, 0.05, 0.12):
            p = delta / alpha
            binary = 1 + p * np.log2(p) + (1 - p) * np.log2(1 - p)
            self.assertAlmostEqual(orbit_rate(np.array([0.0, alpha]), delta), binary, places=9)
        self.assertEqual(orbit_rate(np.array([0.0, alpha]), alpha / 2), 0.0)

    def test_rayleigh_ergodic_closed_form(self):
        for snr in (0.3, 1.0, 10.0, 100.0):
            numeric, _ = integrate.quad(lambda g: np.log2(1 + g * snr) * np.exp(-g), 0, np.inf)
            self.assertAlmostEqual(numeric, np.exp(1 / snr) * exp1(1 / snr) / np.log(2), places=7)


if __name__ == "__main__":
    unittest.main()
