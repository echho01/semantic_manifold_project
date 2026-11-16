"""Finite-alphabet semantic rate-distortion with shared side information.

Costs are fixed by the caller. Rates and the BA certificate use bits per source
symbol. A slope beta has units bits per unit of the supplied distortion.
"""
from dataclasses import dataclass
from itertools import combinations

import numpy as np


@dataclass(frozen=True)
class OperatingPoint:
    beta: float
    rate: float
    distortion: float
    dual: float  # Lower bound on min(I + beta * D) for these alphabets.
    iterations: int

    @property
    def gap(self):
        return max(0.0, self.rate + self.beta * self.distortion - self.dual)


def _mass(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not np.isfinite(values).all() or np.any(values < 0)
            or not np.isclose(values.sum(), 1.0, rtol=0, atol=1e-8)):
        raise ValueError("expected a finite nonnegative probability vector summing to one")
    return values / values.sum()


def finite_ba(mass, cost, beta, tol=1e-7, max_iter=10000):
    """Compute an achievable point and a dual certificate, even at an iteration cap."""
    p = _mass(mass)
    c = np.asarray(cost, dtype=float)
    if (c.ndim != 2 or c.shape[0] != len(p) or c.shape[1] == 0
            or not np.isfinite(c).all() or np.any(c < 0)):
        raise ValueError("expected a finite nonnegative source-by-reproduction cost matrix")
    if (not np.isfinite(beta) or beta < 0 or not np.isfinite(tol) or tol <= 0
            or not isinstance(max_iter, (int, np.integer)) or max_iter < 1):
        raise ValueError("invalid slope, tolerance, or iteration limit")
    if beta == 0:
        return OperatingPoint(0.0, 0.0, float(np.min(p @ c)), 0.0, 0)

    c, p = c[p > 0], p[p > 0]
    row_min = c.min(axis=1)
    kernel = np.exp2(-beta * (c - row_min[:, None]))
    prior = np.full(c.shape[1], 1.0 / c.shape[1])
    # Subtracting the minimum cost of each row keeps a kernel entry of one.
    # Maintain support for columns that can become competitive after an update.
    for iteration in range(1, max_iter + 1):
        z = kernel @ prior
        column_gain = (p / z) @ kernel
        log_gain = np.log2(column_gain.max())
        if log_gain <= tol or iteration == max_iter:
            break
        prior = np.maximum(prior * column_gain, 1e-250)
        prior /= prior.sum()

    transition = (kernel * prior) / z[:, None]
    joint = p[:, None] * transition
    marginal = joint.sum(axis=0)
    log_channel = np.log2(transition, where=transition > 0, out=np.zeros_like(transition))
    log_output = np.log2(marginal, where=marginal > 0, out=np.zeros_like(marginal))
    rate = max(0.0, float(np.sum(joint * (log_channel - log_output))))
    distortion = float(np.sum(joint * c))
    dual = float(beta * (p @ row_min) - p @ np.log2(z) - log_gain)
    return OperatingPoint(float(beta), rate, distortion, dual, iteration)


def conditional_ba(context_mass, source_masses, costs, beta, **kwargs):
    """Shared context at both terminals, one average distortion budget and one slope."""
    weights = _mass(context_mass)
    if len(source_masses) != len(weights) or len(costs) != len(weights):
        raise ValueError("expected one source and one cost matrix per context")
    points = [(w, finite_ba(p, c, beta, **kwargs))
              for w, p, c in zip(weights, source_masses, costs) if w > 0]
    return OperatingPoint(float(beta), sum(w * a.rate for w, a in points),
                          sum(w * a.distortion for w, a in points),
                          sum(w * a.dual for w, a in points),
                          max(a.iterations for _, a in points))


def bounds_at_distortion(points, target):
    """Finite RD bracket from supporting lines and time-sharing of feasible channels.

    An infinite upper bound means the supplied channels do not reach the target.
    Pass only operating points of the same source, criterion and reproduction set.
    """
    points = list(points)
    if not points or not np.isfinite(target) or target < 0:
        raise ValueError("expected operating points and a finite nonnegative target")
    lower = max(0.0, max(a.dual - a.beta * target for a in points))
    upper = min((a.rate for a in points if a.distortion <= target), default=np.inf)
    for a, b in combinations(sorted(points, key=lambda pt: pt.distortion), 2):
        if a.distortion < target < b.distortion:
            frac = (target - a.distortion) / (b.distortion - a.distortion)
            upper = min(upper, (1.0 - frac) * a.rate + frac * b.rate)
    return float(lower), float(upper)


def posterior_kl_bits(source, reproduction):
    """Pairwise directed KL for strictly positive finite task posteriors."""
    p, r = (np.asarray(a, dtype=float) for a in (source, reproduction))
    for a in (p, r):
        if (a.ndim != 2 or min(a.shape) < 1 or not np.isfinite(a).all()
                or np.any(a <= 0) or not np.allclose(a.sum(1), 1, rtol=0, atol=1e-7)):
            raise ValueError("expected rows of strictly positive class probabilities")
    if p.shape[1] != r.shape[1]:
        raise ValueError("source and reproduction must use the same task classes")
    return np.maximum((p * np.log2(p)).sum(1)[:, None] - p @ np.log2(r).T, 0.0)
