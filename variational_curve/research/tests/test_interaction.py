"""Random-instance checks of the mass-absorption, stability and interaction-rank statements.

Only NumPy and SciPy are required:  python tests/test_interaction.py
"""
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from interaction_rank import (  # noqa: E402
    LN2,
    blahut,
    center,
    circle,
    distortion,
    gibbs,
    information_nats,
    interaction_sup,
    pythagoras,
    stability_constant_bits,
)

warnings.filterwarnings("ignore", message=".*encountered in matmul")


def random_problem(rng, n=9, m=6):
    p = rng.random(n) + 0.2
    return p / p.sum(), rng.random((n, m)) * 3.0


class Absorption(unittest.TestCase):
    def test_separable_change_is_absorbed_by_masses(self):
        rng = np.random.default_rng(1)
        p, cost = random_problem(rng)
        mass = rng.random(cost.shape[1]) + 0.1
        a, b, s = rng.normal(size=cost.shape[0]), rng.normal(size=cost.shape[1]), 1.7
        moved = cost + a[:, None] + b[None, :]
        same = gibbs(np.log(mass) + s * b, moved, s)
        np.testing.assert_allclose(same, gibbs(np.log(mass), cost, s), atol=1e-12)

    def test_cross_difference_changes_log_odds(self):
        rng = np.random.default_rng(2)
        p, cost = random_problem(rng)
        bump = np.zeros_like(cost)
        bump[0, 0] = 0.4
        s = 2.0
        residual = center(s * bump)
        self.assertGreater(np.abs(residual).max(), 1e-3)
        # log odds of column 0 against column 1 must be row independent to be a mass change
        odds = -s * (bump[:, 0] - bump[:, 1])
        self.assertGreater(np.ptp(odds), 1e-3)


class Stability(unittest.TestCase):
    def test_ratio_distortion_and_information_bounds(self):
        rng = np.random.default_rng(3)
        for _ in range(300):
            p, cost = random_problem(rng, n=int(rng.integers(3, 12)), m=int(rng.integers(2, 9)))
            n, m = cost.shape
            evaluator = rng.random((n, m)) * 2
            eps, s = rng.random() * 0.6, rng.random() * 3 + 0.1
            a, b = rng.normal(size=n), rng.normal(size=m)
            err = rng.uniform(-eps, eps, size=(n, m))
            moved = cost + a[:, None] + b[None, :] + err
            mass = rng.random(m) + 0.05
            w = gibbs(np.log(mass), cost, s)
            w2 = gibbs(np.log(mass) + s * b, moved, s)
            rho = 2 * s * eps
            ratio = w2 / w
            self.assertLessEqual(ratio.max(), np.exp(rho) * (1 + 1e-10))
            self.assertGreaterEqual(ratio.min(), np.exp(-rho) * (1 - 1e-10))
            d, d2 = np.sum(p[:, None] * w * evaluator), np.sum(p[:, None] * w2 * evaluator)
            self.assertLessEqual(d2, np.exp(rho) * d + 1e-12)
            i, i2 = information_nats(p, w), information_nats(p, w2)
            bound = i + (np.exp(rho) - 1) * (i + 2 / np.e) + min(rho, (np.exp(rho) - 1) ** 2)
            self.assertLessEqual(i2, bound + 1e-12)
            lam = rng.random() * 4
            lag = i / LN2 + lam * d
            lag2 = i2 / LN2 + lam * d2
            self.assertLessEqual(lag2, np.exp(rho) * lag + stability_constant_bits(rho) + 1e-12)

    def test_excess_identity_and_quadratic_bound(self):
        rng = np.random.default_rng(8)
        for _ in range(200):
            p, evaluator = random_problem(rng, n=int(rng.integers(3, 12)), m=int(rng.integers(2, 9)))
            n, m = evaluator.shape
            beta = rng.random() * 3 + 0.2
            pi = blahut(p, evaluator, beta, tol=1e-15, maxiter=200000)
            eps = rng.random() * 1.5
            a, b = rng.normal(size=n), rng.normal(size=m)
            cost = evaluator + (a[:, None] + b[None, :] + rng.uniform(-eps, eps, (n, m))) / beta
            out = pythagoras(p, evaluator, cost, beta, pi)
            # exact identity: excess = sum_i p_i KL(V_i||W*_i) - KL(Vbar||pi*)
            self.assertAlmostEqual(out["excess_bits"], out["identity_bits"], places=7)
            self.assertGreaterEqual(out["excess_bits"], -1e-7)
            self.assertLessEqual(out["excess_bits"], out["weighted_divergence_bits"] + 1e-9)
            self.assertLessEqual(out["weighted_divergence_bits"], out["sup_bound_bits"] + 1e-9)
            # an arbitrary channel also obeys the identity-based upper bound
            v = rng.random((n, m)) + 0.01
            v /= v.sum(1, keepdims=True)
            star = gibbs(np.log(np.maximum(pi, 1e-300)), evaluator, beta)
            gap = (
                information_nats(p, v) + beta * distortion(p, v, evaluator)
                - information_nats(p, star) - beta * distortion(p, star, evaluator)
            )
            kl = float(np.sum(p[:, None] * v * (np.log(v) - np.log(np.maximum(star, 1e-300)))))
            self.assertLessEqual(gap, kl + 1e-7)

    def test_linear_program_never_exceeds_centering(self):
        rng = np.random.default_rng(9)
        for _ in range(20):
            delta = rng.normal(size=(6, 5))
            best = interaction_sup(delta)
            self.assertLessEqual(best, np.abs(center(delta)).max() + 1e-9)
            sep = rng.normal(size=6)[:, None] + rng.normal(size=5)[None, :]
            self.assertAlmostEqual(interaction_sup(delta + sep), best, places=7)
            self.assertAlmostEqual(interaction_sup(sep), 0.0, places=7)

    def test_centering_gives_feasible_separable_part(self):
        rng = np.random.default_rng(4)
        delta = rng.normal(size=(7, 5))
        res = center(delta)
        sep = delta - res
        cross = sep[0, 0] - sep[0, 1] - sep[1, 0] + sep[1, 1]
        self.assertAlmostEqual(cross, 0.0, places=12)
        np.testing.assert_allclose(res.sum(0), 0, atol=1e-12)
        np.testing.assert_allclose(res.sum(1), 0, atol=1e-12)


class InteractionRank(unittest.TestCase):
    def test_bilinear_surface_of_rank_dimension_reproduces_blahut_channel(self):
        rng = np.random.default_rng(5)
        for d in (1, 2, 3):
            n, m = 12, 8
            p = rng.random(n) + 0.3
            p /= p.sum()
            low = rng.normal(size=(n, d)) @ rng.normal(size=(m, d)).T
            cost = low + rng.random(n)[:, None] + rng.random(m)[None, :]
            cost -= cost.min()
            t = center(cost)
            self.assertEqual(np.linalg.matrix_rank(t, tol=1e-9), d)
            beta = 0.8  # nat slope, i.e. lambda*ln2
            pi = np.maximum(blahut(p, cost, beta, tol=1e-14), 1e-300)
            target = gibbs(np.log(pi), cost, beta)
            u, sv, vt = np.linalg.svd(t)
            s = 2.3
            x = -(beta / (2 * s)) * u[:, :d] * sv[:d]
            y = vt[:d].T
            surface = ((x[:, None, :] - y[None, :, :]) ** 2).sum(-1)
            # masses that cancel the separable column term
            col = center_parts(s * surface - beta * cost)[1]
            got = gibbs(np.log(pi) + col, surface, s)
            np.testing.assert_allclose(got, target, atol=1e-9)
            # one dimension fewer leaves a nonzero interaction residual
            if d > 1:
                short = (u[:, : d - 1] * sv[: d - 1]) @ vt[: d - 1]
                self.assertGreater(np.abs(t - short).max(), 1e-6)

    def test_rank_is_necessary(self):
        rng = np.random.default_rng(6)
        n, m, d = 10, 7, 2
        cost = rng.random((n, m))
        t = center(cost)
        self.assertGreater(np.linalg.matrix_rank(t, tol=1e-9), d)
        # any d-dimensional bilinear cost has centred rank at most d
        x, y = rng.normal(size=(n, d)), rng.normal(size=(m, d))
        surface = ((x[:, None, :] - y[None, :, :]) ** 2).sum(-1)
        self.assertLessEqual(np.linalg.matrix_rank(center(surface), tol=1e-9), d)


def center_parts(matrix):
    """Return (row term, column term) of the separable part  a_i + b_j  of a separable matrix."""
    grand = matrix.mean()
    return matrix.mean(1) - grand, matrix.mean(0)


class SingleSurface(unittest.TestCase):
    def test_arc_is_not_of_negative_type_but_chord_is(self):
        n = 96
        x = 2 * np.pi * (np.arange(n) + 0.5) / n
        gap = np.abs(x[:, None] - x[None, :])
        arc = np.minimum(gap, 2 * np.pi - gap)
        self.assertLess(np.linalg.eigvalsh(-center(arc**2)).min(), -1.0)
        self.assertGreater(np.linalg.eigvalsh(-center(2 - 2 * np.cos(arc))).min(), -1e-9)
        self.assertGreater(np.linalg.eigvalsh(-center(arc)).min(), -1e-9)

    def test_symmetric_separable_term_must_vanish(self):
        rng = np.random.default_rng(7)
        n = 6
        a, b = rng.normal(size=n), rng.normal(size=n)
        sep = a[:, None] + b[None, :]
        # a symmetric, zero-diagonal separable matrix is zero
        sym_violation = np.abs(sep - sep.T).max() + np.abs(np.diag(sep)).max()
        self.assertGreater(sym_violation, 1e-6)
        zero = np.full(n, 0.3)[:, None] - np.full(n, 0.3)[None, :]
        np.testing.assert_allclose(zero, 0, atol=1e-15)


if __name__ == "__main__":
    unittest.main()
