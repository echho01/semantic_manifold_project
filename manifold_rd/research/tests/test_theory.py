"""Numpy-only checks of the structural results in the Transactions manuscript.

Run: python tests/test_theory.py
"""
import sys
import unittest
from itertools import permutations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from theory_checks import (LN2, blahut, circle, fisher_speed_and_curvature, information_bits,
                           kl_bits, loss_budget)
from manifold_rd import finite_ba, posterior_kl_bits


class LossBudget(unittest.TestCase):
    def test_identity_on_enumerated_link(self):
        for row in loss_budget(beta=5.0):
            self.assertLess(row["identity_error"], 1e-12)
            self.assertGreaterEqual(row["channel"], -1e-12)
            self.assertGreaterEqual(row["decoder"], -1e-12)

    def test_random_chain(self):
        rng = np.random.default_rng(7)
        for _ in range(20):
            ny, nx, nu, nj, nv = 4, 9, 5, 4, 6
            p_x = rng.dirichlet(np.ones(nx))
            q_x = rng.dirichlet(np.ones(ny), size=nx)
            u_of_x = rng.integers(0, nu, nx)
            u_of_x[:nu] = np.arange(nu)
            w = rng.dirichlet(np.ones(nj), size=nu)
            k = rng.dirichlet(np.ones(nv), size=nj)
            b = rng.dirichlet(np.ones(ny), size=nv)
            joint = np.zeros((nx, nj, nv))
            for x in range(nx):
                joint[x] = p_x[x] * w[u_of_x[x]][:, None] * k
            def cond_entropy(group):          # H(Y | group) with group an index array over (x, j, v)
                total = 0.0
                for key in np.unique(group):
                    mass = np.where(group == key, joint, 0.0)
                    weight = mass.sum()
                    if weight > 0:
                        mean = np.tensordot(mass.sum((1, 2)), q_x, axes=1) / weight
                        total += weight * float(-(mean * np.log2(mean)).sum())
                return total
            xs, js, vs = np.meshgrid(np.arange(nx), np.arange(nj), np.arange(nv), indexing="ij")
            h_x, h_u = cond_entropy(xs), cond_entropy(u_of_x[xs])
            h_j, h_v = cond_entropy(js), cond_entropy(vs)
            decoder = 0.0
            for v in range(nv):
                mass = joint[:, :, v].sum(1)
                mean = mass @ q_x / mass.sum()
                decoder += mass.sum() * float(kl_bits(mean, b[v]))
            direct = float(sum(joint[x].sum(0) @ (-(q_x[x][None, :] * np.log2(b)).sum(1)) for x in range(nx)))
            budget = h_x + (h_u - h_x) + (h_j - h_u) + (h_v - h_j) + decoder
            self.assertAlmostEqual(direct, budget, places=11)
            self.assertGreaterEqual(h_u - h_x, -1e-12)
            self.assertGreaterEqual(h_j - h_u, -1e-12)
            self.assertGreaterEqual(h_v - h_j, -1e-12)


class PriceOfAnObject(unittest.TestCase):
    def test_centroid_defect_on_circle_and_sphere(self):
        rng = np.random.default_rng(3)
        for dim in (2, 3, 5):
            for _ in range(200):
                pts = rng.normal(size=(rng.integers(2, 6), dim))
                pts /= np.linalg.norm(pts, axis=1, keepdims=True)     # unit sphere, reach one
                weights = rng.dirichlet(np.ones(len(pts)))
                mean = weights @ pts
                var = float(weights @ ((pts - mean) ** 2).sum(1))
                defect = 1 - np.linalg.norm(mean)
                self.assertLessEqual(defect, 1 - np.sqrt(max(0.0, 1 - var)) + 1e-12)

    def test_separation_and_rounding_bounds(self):
        amplitude = 0.2
        eta, radius = 1 / 3 - amplitude, np.sqrt(1.5) * amplitude
        gamma = (eta / (2 * radius)) ** 2
        rng = np.random.default_rng(11)
        fine = circle(np.linspace(0, 2 * np.pi, 20001), amplitude)
        for _ in range(200):
            theta = rng.uniform(0, 2 * np.pi, rng.integers(2, 6))
            weights = rng.dirichlet(np.ones(len(theta)))
            q = circle(theta, amplitude)
            m = weights @ q
            d = float(weights @ kl_bits(q, np.repeat(m[None], len(q), 0))) * LN2
            price = float(np.min(kl_bits(np.repeat(m[None], len(fine), 0), fine))) * LN2
            self.assertGreaterEqual(price, gamma * d * d - 1e-9)
            if d < radius ** 2:
                self.assertLessEqual(price, d * d / (eta * radius ** 2) + 1e-9)

    def test_local_curvature_constant(self):
        theta0 = np.array([0.7])
        _, kappa = fisher_speed_and_curvature(theta0)
        pts = circle(np.array([0.7 - 0.02, 0.7 + 0.02]))
        m = pts.mean(0, keepdims=True)
        d = float(np.mean(kl_bits(pts, np.repeat(m, 2, 0)))) * LN2
        fine = circle(0.7 + np.linspace(-0.05, 0.05, 200001))
        price = float(np.min(kl_bits(np.repeat(m, len(fine), 0), fine))) * LN2
        self.assertAlmostEqual(price / d**2, 0.5 * kappa[0] ** 2, delta=2e-3 * kappa[0] ** 2)

    def test_convex_alphabet_has_no_price(self):
        q = circle(2 * np.pi * np.arange(12) / 12)
        p = np.full(12, 1 / 12)
        w = blahut(p, posterior_kl_bits(q, q), 6.0)
        joint = p[:, None] * w
        mass = joint.sum(0)
        centroids = (joint[:, mass > 1e-12] / mass[mass > 1e-12]).T @ q
        hull = np.vstack([q, centroids])                       # enlarging by the centroids removes the price
        a = finite_ba(p, posterior_kl_bits(q, q), 6.0)
        b = finite_ba(p, posterior_kl_bits(q, hull), 6.0)
        self.assertLess(b.rate + 6.0 * b.distortion, a.rate + 6.0 * a.distortion - 1e-4)
        self.assertGreater(information_bits(p, w), 0)


class AnalogAllocation(unittest.TestCase):
    def test_closed_form_against_search(self):
        rng = np.random.default_rng(5)
        for _ in range(20):
            lam = np.sort(rng.uniform(0.2, 3.0, 3))[::-1]
            gain = np.sort(rng.uniform(0.3, 2.0, 3))[::-1]
            power, noise = 2.5, 0.3
            best = np.inf
            for perm in permutations(range(3)):
                h = gain[list(perm)]
                closed = noise * (np.sqrt(lam) / h).sum() ** 2 / power
                for _ in range(300):                          # random feasible allocations never beat it
                    share = rng.dirichlet(np.ones(3)) * power
                    a2 = share / lam
                    self.assertGreaterEqual(noise * (1 / (h * h * a2)).sum(), closed - 1e-9)
                best = min(best, closed)
            aligned = noise * (np.sqrt(lam) / gain).sum() ** 2 / power
            self.assertAlmostEqual(best, aligned, places=12)


if __name__ == "__main__":
    unittest.main()
