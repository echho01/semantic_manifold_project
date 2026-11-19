"""Variational coding surfaces: isometry defect, anisotropic criteria and the two regimes.

Three deterministic computations support the letter.

  circle   A coding surface whose induced metric h = w(theta)^2 dtheta^2 differs from the
           semantic metric g = dtheta^2.  The isotropic Gibbs coder acts on the surface with
           one global scale.  Its rate excess over the Blahut optimum at matched D_g is
           compared with the isometry defect Lambda[h].
  regimes  The matched surface (h = g) with uniform masses on a nonuniform source: the excess
           is a density effect that vanishes at high resolution.
  torus    A two-dimensional source whose task ignores one coordinate, scored by the joint
           cost KL + eta * distance^2.  Raw, dimension-only and fully matched surfaces.

    python experiments/surface_functional.py --output results/surface_functional.json \
        --figure ../figures/surface_functional.png

Rates are in bits, distortions in the units of the evaluation cost (nats for the joint cost).
Grid values are computed operating points of finite problems, not certificates.
"""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from scipy import integrate

from interaction_rank import LN2, TINY, circle, distortion, gibbs, information_nats, matched
from interaction_rank import blahut as blahut_reference

warnings.filterwarnings("ignore", message=".*encountered in matmul")


DUAL_GAPS = []  # dual gap in nats of every Blahut run of this process


def blahut(p, cost, beta, start=None, tol=2e-8, maxiter=20000):
    """Output marginal of the Blahut iteration (slope in nats), stopped by its dual gap.

    The kernel is formed once, so an iteration is two matrix-vector products.  The gap between
    the achieved Lagrangian and its dual certificate is at most ln(max gain); it is recorded
    in DUAL_GAPS whether or not the tolerance was met.
    """
    kernel = np.exp(-beta * (cost - cost.min(1, keepdims=True)))
    count = cost.shape[1]
    # A warm start keeps every column alive, otherwise a mass that decayed at a smaller slope
    # could not return.
    pi = np.full(count, 1.0 / count) if start is None else np.maximum(start, 1e-9 / count)
    pi = pi / pi.sum()
    for _ in range(maxiter):
        gain = (p / (kernel @ pi)) @ kernel
        gap = float(np.log(gain.max()))
        if gap <= tol:
            break
        pi = np.maximum(pi * gain, 1e-250)
        pi /= pi.sum()
    DUAL_GAPS.append(gap)
    return pi


def circular_sq(coordinate, period):
    gap = np.abs(coordinate[:, None] - coordinate[None, :])
    return np.minimum(gap, period - gap) ** 2


def reference_curve(p, evaluator, slopes):
    """Blahut operating points (D_g, bits) of the evaluation cost itself."""
    points, pi = [], None
    for beta in slopes:
        pi = blahut(p, evaluator, beta, pi)
        channel = gibbs(np.log(np.maximum(pi, TINY)), evaluator, beta)
        points.append((distortion(p, channel, evaluator), information_nats(p, channel) / LN2))
    return points


def surface_curve(p, coding, evaluator, scales, optimize_mass):
    """Operating points of the isotropic Gibbs coder on a surface, one global scale each."""
    points, pi = [], None
    for scale in scales:
        if optimize_mass:
            pi = blahut(p, coding, scale, pi)
            log_mass = np.log(np.maximum(pi, TINY))
        else:
            log_mass = np.zeros(coding.shape[1])
        channel = gibbs(log_mass, coding, scale)
        points.append((distortion(p, channel, evaluator), information_nats(p, channel) / LN2))
    return points


def log_matched(points, target):
    """Rate at a target distortion, interpolating linearly in ln D."""
    ordered = np.array(sorted(points))
    if not ordered[0, 0] <= target <= ordered[-1, 0]:
        return None
    return float(np.interp(np.log(target), np.log(ordered[:, 0]), ordered[:, 1]))


def gaps(family, reference, targets):
    out = {}
    for target in targets:
        a, b = log_matched(family, target), log_matched(reference, target)
        out[f"{target:g}"] = None if a is None or b is None else a - b
    return out


def defect_1d(weights, stretch_sq):
    """Isometry defect in nats for a curve: h = stretch_sq * g, source weights on a grid."""
    return 0.5 * float(np.log(weights @ (1.0 / stretch_sq)) + weights @ np.log(stretch_sq))


def circle_study(n, tilt, targets):
    theta = 2 * np.pi * (np.arange(n) + 0.5) / n
    evaluator = circular_sq(theta, 2 * np.pi)
    stretch = 1 + tilt * np.cos(theta)                       # sqrt(h/g)
    coding = circular_sq(theta + tilt * np.sin(theta), 2 * np.pi)
    slopes = np.geomspace(0.3, 600, 44)
    out = {"states": n, "stretch": f"1 + {tilt} cos(theta)"}
    for name, modulation in (("uniform_source", 0.0), ("nonuniform_source", 0.85)):
        p = (1 + modulation * np.cos(theta)) / n
        reference = reference_curve(p, evaluator, slopes)
        entry = {
            "defect_bits": defect_1d(p, stretch**2) / LN2,
            "mismatched_uniform_mass": gaps(surface_curve(p, coding, evaluator, slopes, False), reference, targets),
            "mismatched_optimized_mass": gaps(surface_curve(p, coding, evaluator, slopes, True), reference, targets),
            "matched_uniform_mass": gaps(surface_curve(p, evaluator, evaluator, slopes, False), reference, targets),
        }
        out[name] = entry
    closed = 0.5 * np.log((1 - tilt**2) ** -1.5) + np.log((1 + np.sqrt(1 - tilt**2)) / 2)
    out["uniform_source"]["defect_closed_form_bits"] = float(closed / LN2)
    return out


def circle_gibbs(scale):
    """Exact (mean squared arc, information in nats) of exp(-scale x^2) on the unit circle."""
    cut = min(np.pi, 12.0 / np.sqrt(scale))
    z = 2 * integrate.quad(lambda x: np.exp(-scale * x * x), 0, cut, epsabs=0, epsrel=1e-13)[0]
    second = 2 * integrate.quad(lambda x: x * x * np.exp(-scale * x * x), 0, cut, epsabs=0, epsrel=1e-13)[0] / z
    return second, np.log(2 * np.pi) - np.log(z) - scale * second


def torus_study(n, amplitude, eta, targets):
    theta = 2 * np.pi * (np.arange(n) + 0.5) / n
    phases = 2 * np.pi * np.arange(3) / 3
    q = 1 / 3 + amplitude * np.cos(theta[:, None] - phases)
    dq = -amplitude * np.sin(theta[:, None] - phases)
    fisher = (dq * dq / q).sum(1)                             # squared Fisher speed
    kl = (q * np.log(q)).sum(1)[:, None] - q @ np.log(q).T    # nats
    arc2 = circular_sq(theta, 2 * np.pi)
    evaluator = np.maximum(kl, 0.0) + eta * arc2              # theta part of the joint cost
    p = np.full(n, 1 / n)
    sensitive = 0.5 * fisher + eta                            # g_eta along theta in the raw chart
    mean_sensitive = float(sensitive.mean())
    # Arclength of the semantic metric along theta, on the same grid (midpoint rule).
    speed = np.sqrt(sensitive)
    length = float(speed.sum() * 2 * np.pi / n)
    arclength = np.cumsum(speed) * 2 * np.pi / n - 0.5 * speed * 2 * np.pi / n
    surfaces = {
        "raw": (arc2, 1.0),
        "dimension_only": (arc2, eta / mean_sensitive),
        "matched": (circular_sq(arclength, length), eta),
    }
    slopes = np.geomspace(2, 6000, 44)                        # nats of information per nat of distortion
    reference, pi = [], None
    for beta in slopes:
        pi = blahut(p, evaluator, beta, pi)
        channel = gibbs(np.log(np.maximum(pi, TINY)), evaluator, beta)
        second, info = circle_gibbs(beta * eta)
        reference.append((distortion(p, channel, evaluator) + eta * second,
                          (information_nats(p, channel) + info) / LN2))
    # Predicted defects, in nats: eigenvalues of h^{-1} g are (sensitive, eta / eps^2).
    mean_log = float(np.log(sensitive).mean())
    predicted = {
        "raw": np.log((mean_sensitive + eta) / 2) - 0.5 * mean_log - 0.5 * np.log(eta),
        "dimension_only": 0.5 * np.log(mean_sensitive) - 0.5 * mean_log,
        "matched": 0.0,
    }
    out = {"states": n, "eta": eta, "mean_sensitive_eigenvalue": mean_sensitive,
           "optimal_null_stretch_sq": eta / mean_sensitive, "surfaces": {}}
    curves = {"reference": reference}
    for name, (coding, null_stretch_sq) in surfaces.items():
        scales = np.geomspace(0.5, 3e5, 60) if name != "matched" else slopes
        points = []
        for scale in scales:
            channel = gibbs(np.zeros(n), coding, scale)
            second, info = circle_gibbs(scale * null_stretch_sq)
            points.append((distortion(p, channel, evaluator) + eta * second,
                           (information_nats(p, channel) + info) / LN2))
        curves[name] = points
        out["surfaces"][name] = {"defect_bits": float(predicted[name] / LN2),
                                 "gap_bits": gaps(points, reference, targets)}
    return out, curves


def circle_curves(n, tilt):
    """Gap curves for the figure: (D, gap in bits) on a common distortion grid."""
    theta = 2 * np.pi * (np.arange(n) + 0.5) / n
    evaluator = circular_sq(theta, 2 * np.pi)
    coding = circular_sq(theta + tilt * np.sin(theta), 2 * np.pi)
    slopes = np.geomspace(0.3, 600, 44)
    grid = np.geomspace(4e-3, 1.0, 40)
    out = {"distortion": grid.tolist()}
    uniform = np.full(n, 1 / n)
    tilted = (1 + 0.85 * np.cos(theta)) / n
    ref_u, ref_t = reference_curve(uniform, evaluator, slopes), reference_curve(tilted, evaluator, slopes)
    for name, p, cost, ref, opt in (("mismatched", uniform, coding, ref_u, False),
                                    ("mismatched_optimized", uniform, coding, ref_u, True),
                                    ("matched_uniform_mass", tilted, evaluator, ref_t, False)):
        family = surface_curve(p, cost, evaluator, slopes, opt)
        out[name] = [None if (a := log_matched(family, d)) is None or (b := log_matched(ref, d)) is None
                     else a - b for d in grid]
    return out


def plot(circle_gap, circle_result, torus_result, torus_curves, path):
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(3.5, 1.9), layout="constrained")
    d = np.array(circle_gap["distortion"])
    for key, color, label in (("mismatched", "#124971", r"$h\neq g$, uniform mass"),
                              ("mismatched_optimized", "#6a8f3c", r"$h\neq g$, optimized mass"),
                              ("matched_uniform_mass", "#a32437", r"$h=g$, uniform mass")):
        y = np.array([np.nan if v is None else v for v in circle_gap[key]])
        ax.semilogx(d, y, color=color, lw=1.1, label=label)
    ax.axhline(circle_result["uniform_source"]["defect_bits"], color=".4", ls=":", lw=.8)
    ax.set(xlabel=r"$D_g$", ylabel="Excess rate (bits)", title="(a) Circle", ylim=(-0.01, 0.36))
    ax.legend(frameon=False, fontsize=4.6, loc="upper left")
    ref = torus_curves["reference"]
    grid = np.geomspace(4e-4, 0.2, 40)
    for key, color, label in (("raw", "#124971", "raw"), ("dimension_only", "#6a8f3c", "dimension only"),
                              ("matched", "#a32437", "matched")):
        y = [np.nan if (a := log_matched(torus_curves[key], t)) is None or (b := log_matched(ref, t)) is None
             else a - b for t in grid]
        bx.semilogx(grid, y, color=color, lw=1.1, label=label)
        bx.axhline(torus_result["surfaces"][key]["defect_bits"], color=color, ls=":", lw=.7)
    bx.set(xlabel=r"$D$ (nats)", title="(b) Torus", ylim=(-0.05, 1.3))
    bx.legend(frameon=False, fontsize=4.6, loc="upper left")
    for axis in (ax, bx):
        axis.grid(alpha=.2, which="both")
        axis.title.set_fontsize(7)
        axis.tick_params(labelsize=5.5)
        axis.xaxis.label.set_fontsize(6.5)
        axis.yaxis.label.set_fontsize(6.5)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=400)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/surface_functional.json"))
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--states", type=int, default=720)
    args = parser.parse_args()

    # The same finite solver reproduces the reference of the 96-state study.
    p96, arc96 = circle()
    betas = np.geomspace(0.25, 32, 25)
    early = []
    for beta in betas:
        channel = gibbs(np.log(np.maximum(blahut_reference(p96, arc96**2, beta), TINY)), arc96**2, beta)
        early.append((distortion(p96, channel, arc96**2), information_nats(p96, channel) / LN2))
    mine = matched(reference_curve(p96, arc96**2, betas), 0.30)
    early = matched(early, 0.30)
    if abs(mine - early) > 2e-6:
        raise ArithmeticError(f"fast Blahut iteration disagrees with the reference solver: {mine} vs {early}")
    if abs(early - 1.183191) > 2e-6:
        raise ArithmeticError(f"Blahut reference of the 96-state circle changed: {early}")

    circle_targets = (1.0, 0.3, 0.1, 0.03, 0.01, 0.005)
    torus_targets = (0.1, 0.03, 0.01, 0.003, 0.001, 0.0005)
    circle_result = circle_study(args.states, 0.5, circle_targets)
    torus_result, torus_curves = torus_study(args.states, 0.2, 0.01, torus_targets)
    result = {"blahut_96_state_bits_at_0.30": early, "max_blahut_dual_gap_nats": max(DUAL_GAPS),
              "circle": circle_result, "torus": torus_result}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    if args.figure:
        plot(circle_curves(args.states, 0.5), circle_result, torus_result, torus_curves, args.figure)
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
