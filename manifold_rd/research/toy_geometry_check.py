"""Reproduce the posterior-circle example and finite-alphabet RD bounds.

Run: python toy_geometry_check.py [--grid 24] [--output DIRECTORY]
"""
import argparse
import csv
from dataclasses import asdict
from math import erfc
from pathlib import Path

import numpy as np

from manifold_rd import bounds_at_distortion, conditional_ba, finite_ba, posterior_kl_bits


def circle(n, amplitude=0.2):
    angle = 2 * np.pi * np.arange(n) / n
    phases = 2 * np.pi * np.arange(3) / 3
    return 1 / 3 + amplitude * np.cos(angle[:, None] - phases)


def zero_rate():
    rows = []
    for amplitude in (0.05, 0.10, 0.20, 0.25, 0.30):
        q = circle(65536, amplitude)
        center = np.full((1, 3), 1 / 3)
        free = posterior_kl_bits(q, center).mean()
        extra = -np.log2(1 - 27 * amplitude**2 * (1 - amplitude) / 4) / 3
        # Angle zero maximizes the product of the three reproduction probabilities.
        restricted = posterior_kl_bits(q, q[:1]).mean()
        np.testing.assert_allclose(restricted - free, extra, atol=5e-13, rtol=0)
        rows.append(dict(amplitude=amplitude, D_pred_bits=free,
                         gap_bits=extra, D_circle_bits=restricted))
    return rows


def binary_check():
    p = [0.5, 0.5]
    cost = [[0, 1], [1, 0]]
    for beta in (0, 0.5, 2, 8):
        a = finite_ba(p, cost, beta)
        d = 1 / (1 + 2**beta)
        expected_rate = 1 + d * np.log2(d) + (1-d) * np.log2(1-d)
        np.testing.assert_allclose([a.distortion, a.rate], [d, expected_rate], atol=1e-12)
        np.testing.assert_allclose(a.gap, 0, atol=1e-12)
    # Context reveals the first of two independent fair bits. Hamming loss
    # measures the two-bit state, using the same four reproductions in each case.
    states = np.arange(4)
    c = ((states[:, None] // 2 != states // 2).astype(float)
         + (states[:, None] % 2 != states % 2))
    contexts = ([0.5, 0.5, 0, 0], [0, 0, 0.5, 0.5])
    known = conditional_ba(p, contexts, [c, c], 2)
    unknown = finite_ba(np.full(4, 0.25), c, 2)
    np.testing.assert_allclose([known.distortion, unknown.distortion], [0.2, 0.4], atol=1e-7)
    one_bit_rate = finite_ba(p, cost, 2).rate
    np.testing.assert_allclose([known.rate, unknown.rate], [one_bit_rate, 2 * one_bit_rate], atol=1e-7)


def rd_grid(n):
    q = circle(n)
    p = np.full(n, 1 / n)
    center = np.full((1, 3), 1 / 3)
    disk = np.vstack([center] + [center + scale * (q - center) for scale in (0.25, 0.5, 0.75, 1)])
    angle = 2 * np.pi * np.arange(n) / n
    step = np.abs(angle[:, None] - angle)
    arc_squared = np.minimum(step, 2*np.pi - step)**2
    # With K shared, half the contexts favour each of two opposite semicircles.
    # Their equally weighted mixture is the same uniform marginal p.
    tilted = [(1 + sign * 0.8 * np.cos(angle)) / n for sign in (-1, 1)]
    circle_cost = posterior_kl_bits(q, q)
    disk_cost = posterior_kl_bits(q, disk)
    slopes = (0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 16, 24, 32, 48, 64, 96)
    rows, brackets = [], []
    for name, weights, masses, costs, targets in [
        ("squared_unit_circle_arc", [1], [p], [arc_squared], (0.1, 0.5)),
        ("posterior_KL_circle", [1], [p], [circle_cost], (0.05, 0.10)),
        ("posterior_KL_disk_grid", [1], [p], [disk_cost], (0.05, 0.10)),
        ("posterior_KL_circle_shared_K", [0.5, 0.5], tilted, [circle_cost]*2, (0.05, 0.10)),
    ]:
        points = []
        for beta in slopes:
            a = conditional_ba(weights, masses, costs, beta, tol=2e-6)
            points.append(a)
            rows.append(dict(grid=n, criterion=name, **asdict(a), gap_bits=a.gap))
        for target in targets:
            lower, upper = bounds_at_distortion(points, target)
            if not (np.isfinite(upper) and lower <= upper + 1e-10):
                raise ArithmeticError("invalid finite-alphabet RD bracket")
            brackets.append(dict(grid=n, criterion=name, target_distortion=target,
                                 lower_rate_bits=lower, upper_rate_bits=upper))
    return rows, brackets


def gaussian_link(n):
    """Uniform grid sent as (cos theta, sin theta); receiver knows scalar gain.

    Two real AWGN uses, unit noise variance per use, total symbol energy 4.
    Positive gains use maximum likelihood. At zero gain use the best constant
    task reproduction. Integrate angular decision bins of a shifted 2-D Gaussian.
    """
    q = circle(n)
    task_cost = posterior_kl_bits(q, q)
    angle = 2 * np.pi * np.arange(n) / n
    arc_cost = np.minimum(angle, 2 * np.pi - angle)**2
    nodes, weights = np.polynomial.legendre.leggauss(48)
    phi = angle[:, None] + (np.pi / n) * nodes
    offsets = (np.arange(n)[None, :] - np.arange(n)[:, None]) % n
    rows = []
    for gain in (0.0, 0.25, 0.5, 1.0, 2.0):
        energy = 4.0
        rho = gain * np.sqrt(energy)
        if gain == 0:
            channel = np.zeros((n, n))
            channel[:, np.argmin(task_cost.mean(axis=0))] = 1
        else:
            a = rho * np.cos(phi)
            normal_cdf = np.array([0.5 * erfc(-v/np.sqrt(2)) for v in a.flat]).reshape(a.shape)
            density = (np.exp(-rho**2 / 2) / (2 * np.pi)
                       + a * np.exp((a*a - rho*rho) / 2) * normal_cdf / np.sqrt(2*np.pi))
            probs = (np.pi / n) * (density @ weights)
            np.testing.assert_allclose(probs.sum(), 1.0, atol=1e-12, rtol=0)
            if np.min(probs) < 0:
                raise ArithmeticError("negative angular bin probability")
            channel = probs[offsets]
        marginal = channel.mean(axis=0)
        rate = float(np.sum(channel * np.log2(channel, where=channel > 0, out=np.zeros_like(channel))) / n
                     - np.sum(marginal * np.log2(marginal, where=marginal > 0, out=np.zeros_like(marginal))))
        capacity = float(np.log2(1 + gain**2 * energy / 2))
        if rate > capacity + 1e-10:
            raise ArithmeticError("decoded mutual information exceeds the two-use capacity")
        rows.append(dict(grid=n, gain=gain, total_energy=energy, real_uses=2,
                         capacity_bits_per_state=capacity, decoded_MI_bits=max(0.0, rate),
                         D_arc_squared=float(np.sum(channel * arc_cost[offsets]) / n),
                         D_task_KL_bits=float(np.sum(channel * task_cost) / n)))
    return rows


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def plot_posterior_circle(output):
    """Reproduce the zero-rate figure directly from the circle source and its KL centroid."""
    import matplotlib.pyplot as plt

    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    points = circle(360) @ vertices
    center = vertices.mean(0)
    fig, (plane, ax) = plt.subplots(1, 2, figsize=(6.9, 2.65), layout="constrained",
                                    gridspec_kw={"width_ratios": [1, 1.25]})
    plane.plot(*np.vstack((vertices, vertices[0])).T, color=".55", lw=1)
    plane.fill(*points.T, color="#e2ecf2")
    plane.plot(*np.vstack((points, points[0])).T, color="#124971", lw=1.7)
    plane.plot([center[0], points[0, 0]], [center[1], points[0, 1]], "--", color=".4", lw=1)
    plane.scatter(*center, color="#a32437", s=17, zorder=3)
    plane.scatter(*points[0], color="#124971", s=17, zorder=3)
    plane.annotate("Mean posterior", center, (.57, .54), fontsize=8,
                   arrowprops={"arrowstyle": "-", "lw": .7})
    plane.annotate("Best object", points[0], (-.02, .39), fontsize=8,
                   arrowprops={"arrowstyle": "-", "lw": .7})
    for (x, y), label, offset in zip(vertices, ("Class 1", "Class 2", "Class 3"), (-10, -10, 4)):
        plane.annotate(label, (x, y), xytext=(0, offset), textcoords="offset points",
                       fontsize=7, ha="center")
    plane.set(xlim=(-.12, 1.12), ylim=(-.11, 1), aspect="equal", title="(a) Posterior circle")
    plane.axis("off")
    amplitudes = np.linspace(0, .31, 112)
    c = np.full((1, 3), 1/3)
    free = np.array([posterior_kl_bits(circle(4096, a), c).mean() for a in amplitudes])
    penalty = -np.log2(1 - 27 * amplitudes**2 * (1 - amplitudes) / 4) / 3
    ax.fill_between(amplitudes, free, free + penalty, color="#e2ecf2")
    ax.plot(amplitudes, free + penalty, color="#124971", label="Object")
    ax.plot(amplitudes, free, color="#a32437", label="Free predictor")
    ax.set(xlabel="Amplitude $a$", ylabel="Zero-rate KL (bits)", xlim=(0, .31), ylim=(0, None),
           title="(b) Receiver constraint")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=.2)
    for axis in (plane, ax):
        axis.title.set_fontsize(9)
        axis.tick_params(labelsize=7)
    ax.xaxis.label.set_fontsize(8)
    ax.yaxis.label.set_fontsize(8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=24)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--figure", type=Path, help="write the manuscript's zero-rate figure")
    args = parser.parse_args()
    if args.grid < 6 or args.grid % 3:
        parser.error("--grid must be a multiple of 3 and at least 6")
    binary_check()
    zero, (grid, bounds) = zero_rate(), rd_grid(args.grid)
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "toy_zero_rate.csv", zero)
    write_csv(args.output / "finite_rd.csv", grid)
    write_csv(args.output / "finite_rd_bounds.csv", bounds)
    write_csv(args.output / "gaussian_link.csv", gaussian_link(args.grid))
    if args.figure:
        plot_posterior_circle(args.figure)
    print(f"Binary checks passed. Wrote {len(grid)} operating points and {len(bounds)} bounds to {args.output}")


if __name__ == "__main__":
    main()
