"""Numpy-only checks of the variational surface results.

Run: python tests/test_variational.py
"""
import sys
import unittest
from math import gamma, pi
from pathlib import Path

import numpy as np

np.seterr(divide="ignore", over="ignore", invalid="ignore")  # spurious Accelerate matmul flags
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))


def spd(rng, m, count):
    a = rng.normal(size=(count, m, m))
    return a @ np.swapaxes(a, 1, 2) + 0.3 * np.eye(m)


def defect(weights, h, g):
    """Isometry defect in nats for metric fields sampled at weighted points."""
    m = h.shape[-1]
    ratio = np.linalg.solve(h, g)
    trace = np.trace(ratio, axis1=1, axis2=2)
    logdet = np.linalg.slogdet(h)[1] - np.linalg.slogdet(g)[1]
    return 0.5 * m * np.log(weights @ trace / m) + 0.5 * weights @ logdet


def ball(k):
    return pi ** (k / 2) / gamma(k / 2 + 1)


class IsometryDefect(unittest.TestCase):
    def test_nonnegative_and_zero_only_for_scaled_metric(self):
        rng = np.random.default_rng(1)
        for m in (1, 2, 4):
            g = spd(rng, m, 50)
            weights = rng.dirichlet(np.ones(50))
            for _ in range(20):
                self.assertGreater(defect(weights, spd(rng, m, 50), g), 1e-6)
            self.assertAlmostEqual(defect(weights, 3.7 * g, g), 0.0, places=12)

    def test_sandwich_bound(self):
        rng = np.random.default_rng(2)
        m = 3
        for _ in range(50):
            g = spd(rng, m, 40)
            weights = rng.dirichlet(np.ones(40))
            a, b = 0.4, 2.5
            root = np.linalg.cholesky(g)
            basis = np.linalg.qr(rng.normal(size=(40, m, m)))[0]
            eig = rng.uniform(a, b, size=(40, m))
            h = root @ (basis * eig[:, None, :]) @ np.swapaxes(basis, 1, 2) @ np.swapaxes(root, 1, 2)
            self.assertLessEqual(defect(weights, h, g), 0.5 * m * np.log(b / a) + 1e-12)

    def test_first_variation_by_finite_differences(self):
        rng = np.random.default_rng(3)
        m, count = 3, 30
        g, h = spd(rng, m, count), spd(rng, m, count)
        weights = rng.dirichlet(np.ones(count))
        dh = rng.normal(size=(count, m, m))
        dh = dh + np.swapaxes(dh, 1, 2)
        ratio = np.linalg.solve(h, g)
        mean_trace = weights @ np.trace(ratio, axis1=1, axis2=2)
        inner = np.linalg.solve(h, dh) @ (np.eye(m) - m * ratio / mean_trace)
        formula = 0.5 * weights @ np.trace(inner, axis1=1, axis2=2)
        step = 1e-6
        numeric = (defect(weights, h + step * dh, g) - defect(weights, h - step * dh, g)) / (2 * step)
        self.assertAlmostEqual(formula, numeric, places=8)

    def test_optimal_null_stretch(self):
        rng = np.random.default_rng(4)
        m, r, eta = 5, 2, 0.01
        sensitive = rng.uniform(0.05, 1.0, size=(200, r))
        weights = rng.dirichlet(np.ones(200))
        mean_sensitive = weights @ sensitive.sum(1) / r

        def value(null):                       # null = eta / eps^2, eigenvalue on task-null directions
            trace = weights @ sensitive.sum(1) + (m - r) * null
            return 0.5 * m * np.log(trace / m) - 0.5 * weights @ np.log(sensitive).sum(1) - 0.5 * (m - r) * np.log(null)

        grid = np.geomspace(1e-3, 10, 20001)
        best = grid[np.argmin([value(x) for x in grid])]
        self.assertAlmostEqual(best / mean_sensitive, 1.0, delta=2e-3)
        # A surface that keeps unit stretch on the null directions pays this much more.
        penalty = (0.5 * (m - r) * np.log(mean_sensitive / eta)
                   + 0.5 * m * np.log((r * mean_sensitive + (m - r) * eta) / (m * mean_sensitive)))
        self.assertAlmostEqual(value(eta) - value(mean_sensitive), penalty, places=10)
        self.assertGreater(penalty, 0)

    def test_circle_closed_form(self):
        b = 0.5
        theta = 2 * np.pi * (np.arange(4000) + 0.5) / 4000
        stretch_sq = (1 + b * np.cos(theta)) ** 2
        weights = np.full(4000, 1 / 4000)
        numeric = 0.5 * np.log(weights @ (1 / stretch_sq)) + 0.5 * weights @ np.log(stretch_sq)
        closed = 0.5 * np.log((1 - b * b) ** -1.5) + np.log((1 + np.sqrt(1 - b * b)) / 2)
        self.assertAlmostEqual(numeric, closed, places=10)


class AreaReach(unittest.TestCase):
    @staticmethod
    def volume_cap(k, m, rho, tau):
        return 2**k * ball(k) * (rho + tau / 2) ** k / (ball(k - m) * tau ** (k - m))

    def test_tube_packing_on_explicit_surfaces(self):
        rho = 1.3
        self.assertLessEqual(2 * pi * rho, self.volume_cap(2, 1, rho, rho))                  # circle in the plane
        self.assertLessEqual(4 * pi * rho**2, self.volume_cap(3, 2, rho, rho))               # sphere in R^3
        a = 0.8
        self.assertLessEqual(4 * pi**2 * a**2, self.volume_cap(4, 2, a * np.sqrt(2), a))     # Clifford torus in R^4
        for major, minor in ((2.0, 1.0), (3.0, 0.5)):                                       # ellipse, reach b^2 / a
            t = 2 * np.pi * (np.arange(200000) + 0.5) / 200000
            length = np.sqrt((major * np.sin(t)) ** 2 + (minor * np.cos(t)) ** 2).sum() * 2 * np.pi / 200000
            self.assertLessEqual(length, self.volume_cap(2, 1, major, minor**2 / major))

    def test_tube_volume_lower_bound_on_circle(self):
        # Exact tube area of a circle of radius R is 2 * (2 pi R) * t; the bound uses (1 - t / tau)^m.
        radius, t = 1.0, 0.5
        self.assertGreaterEqual(4 * pi * radius * t, 2 * pi * radius * ball(1) * t * (1 - t / radius))

    def test_weak_noise_distortion_and_lower_bound(self):
        rng = np.random.default_rng(6)
        # Ellipse F(theta) = (a cos, b sin); semantic metric g = dtheta^2, source p on theta.
        a, b, sigma = 2.0, 1.0, 0.02
        theta = 2 * np.pi * (np.arange(20000) + 0.5) / 20000
        surface = np.c_[a * np.cos(theta), b * np.sin(theta)]
        h = (a * np.sin(theta)) ** 2 + (b * np.cos(theta)) ** 2
        p = 1 + 0.5 * np.cos(theta)
        p /= p.sum()
        predicted = sigma**2 * p @ (1 / h)
        sample = rng.choice(len(theta), size=4000, p=p)
        noisy = surface[sample] + sigma * rng.normal(size=(4000, 2))
        errors = []
        for index, point in zip(sample, noisy):
            window = (index + np.arange(-400, 401)) % len(theta)
            nearest = window[np.argmin(((surface[window] - point) ** 2).sum(1))]
            gap = abs(theta[nearest] - theta[index])
            errors.append(min(gap, 2 * np.pi - gap) ** 2)
        measured = np.mean(errors)
        conditional = sigma**2 * np.mean(1 / h[sample])
        self.assertAlmostEqual(measured / conditional, 1.0, delta=0.08)
        density = p * len(theta) / (2 * np.pi)                    # density relative to dtheta
        entropy_g = -(p * np.log(density)).sum()
        length = np.sqrt(h).sum() * 2 * np.pi / len(theta)
        self.assertGreaterEqual(predicted, sigma**2 * np.exp(2 * entropy_g) / length**2 - 1e-12)

    def test_companding_bound_and_its_minimizer(self):
        rng = np.random.default_rng(8)
        for m in (1, 2, 3):
            alpha = m / (m + 2)
            p = rng.dirichlet(np.ones(400)) * 400            # density on a space of total g-volume one
            cell = np.full(400, 1 / 400)
            norm = (cell @ p**alpha) ** (1 / alpha)          # ||p||_{m/(m+2)}
            entropy = -(cell * p) @ np.log(p)
            self.assertGreaterEqual(norm, np.exp(2 * entropy / m) - 1e-12)
            for _ in range(30):                              # conformal factors h = c g
                c = rng.uniform(0.2, 5.0, 400)
                volume = cell @ c ** (m / 2)
                self.assertGreaterEqual((cell * p) @ (m / c), m * norm * volume ** (-2 / m) - 1e-10)
            c = p ** (2 / (m + 2))                           # the minimizer
            volume = cell @ c ** (m / 2)
            self.assertAlmostEqual((cell * p) @ (m / c), m * norm * volume ** (-2 / m), places=10)

    def test_bound_is_tight_for_arclength_uniform_source(self):
        theta = 2 * np.pi * (np.arange(20000) + 0.5) / 20000
        h = (2.0 * np.sin(theta)) ** 2 + np.cos(theta) ** 2
        g = h.copy()                                             # isometric surface
        p = np.sqrt(h) / np.sqrt(h).sum()                        # uniform in surface length
        density = p * len(theta) / (2 * np.pi) / np.sqrt(g)      # relative to semantic length
        entropy_g = -(p * np.log(density)).sum()
        length = np.sqrt(h).sum() * 2 * np.pi / len(theta)
        self.assertAlmostEqual(p @ (g / h), np.exp(2 * entropy_g) / length**2, places=10)

    def test_scalar_curvature_caps_reach(self):
        # Round sphere of radius R in R^{m+1}: Scal = m (m - 1) / R^2 and reach R, so the cap is attained.
        for m, radius in ((2, 0.7), (3, 1.9), (5, 1.0)):
            scal = m * (m - 1) / radius**2
            self.assertAlmostEqual(np.sqrt(m * (m - 1) / scal), radius, places=12)
        # Ellipsoid of revolution x^2/a^2 + (y^2 + z^2)/b^2 = 1 with a > b: at the poles of the long axis
        # both principal curvatures equal a / b^2, so Scal = 2 a^2 / b^4 and the cap sqrt(2 / Scal) = b^2 / a
        # is the true reach.
        a_axis, b_axis = 2.0, 1.0
        self.assertAlmostEqual(np.sqrt(2 / (2 * a_axis**2 / b_axis**4)), b_axis**2 / a_axis, places=12)


class AnalogMapping(unittest.TestCase):
    def test_companded_circle_meets_the_area_bound(self):
        from analog_mapping import exact_distortion
        tilt, snr = 0.85, 1e3
        p = lambda t: (1 + tilt * np.cos(t)) / (2 * np.pi)
        theta = (np.arange(100000) + 0.5) * 2 * np.pi / 100000
        norm = (np.sum(p(theta) ** (1 / 3)) * 2 * np.pi / 100000) ** 3
        bound = norm / (2 * np.pi) ** 2 / snr
        isometric, companded = (exact_distortion(p, e, snr, n=2000) for e in (0.0, 1 / 3))
        self.assertGreaterEqual(companded, bound * (1 - 1e-3))       # the bound holds
        self.assertAlmostEqual(companded / bound, 1.0, delta=0.01)   # and is met by the p^(1/3) law
        self.assertAlmostEqual(isometric * snr, 1.0, delta=0.01)     # isometric surface: sigma^2 / rho^2
        self.assertLess(companded, isometric)


class TwoRegimes(unittest.TestCase):
    def test_defect_limit_and_density_regime(self):
        from surface_functional import LN2, circle_study
        out = circle_study(240, 0.5, (0.3, 0.03))
        uniform = out["uniform_source"]
        self.assertAlmostEqual(uniform["defect_bits"], uniform["defect_closed_form_bits"], places=6)
        for family in ("mismatched_uniform_mass", "mismatched_optimized_mass"):
            self.assertAlmostEqual(uniform[family]["0.03"], uniform["defect_bits"], delta=0.02)
        nonuniform = out["nonuniform_source"]
        self.assertGreater(nonuniform["matched_uniform_mass"]["0.3"], 0.05)
        self.assertLess(nonuniform["matched_uniform_mass"]["0.03"], 0.02)
        self.assertGreater(LN2, 0)


if __name__ == "__main__":
    unittest.main()
