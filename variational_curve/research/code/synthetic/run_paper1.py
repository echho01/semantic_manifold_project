#!/usr/bin/env python3
"""Reproduce Paper 1's numerical tables and figures on a CPU.

All rates are in nats. Reconstruction points lie on the stated manifold.
The sphere curves are continuous Gibbs integrals. The surface curves use
quadrature in the meridian and analytic integration in the azimuth, with
a certified Blahut--Arimoto Lagrangian gap. No learned-model data are used.

Run from the repository root:
    python code/synthetic/run_paper1.py --output code/results/paper1
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
from pathlib import Path
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/semantic-manifold-matplotlib")
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import (gammaln, i0e, i1e, roots_genlaguerre,
                           roots_hermitenorm, roots_legendre)
from scipy.stats import ncx2, norm
from threadpoolctl import threadpool_limits


def extrapolate(x, y, degree=2):
    """Return the zero-scale intercept, with scaling for conditioning."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    degree = min(degree, len(x) - 1)
    scale = np.max(x)
    coeff = np.polynomial.polynomial.polyfit(x / scale, y, degree)
    return float(coeff[0])


def sphere_rd(d, mu, criterion="geodesic"):
    """Uniform S^d(1), exact continuous RD at the specified temperature."""
    if criterion not in ("geodesic", "chordal") or d < 1 or mu <= 0:
        raise ValueError("invalid sphere RD parameters")
    # u=sqrt(mu)*theta keeps the integration window resolved as mu grows.
    hi = min(np.pi * np.sqrt(mu), 14.0)
    def kernel(u, moment):
        theta = u / np.sqrt(mu)
        scaled_dist = u * u if criterion == "geodesic" else 4 * mu * np.sin(theta / 2) ** 2
        jac = u ** (d - 1) * np.sinc(theta / np.pi) ** (d - 1)
        return np.exp(-scaled_dist) * jac * (scaled_dist if moment else 1.0)
    z = quad(lambda u: kernel(u, False), 0, hi, epsabs=1e-13, epsrel=1e-12)[0]
    m = quad(lambda u: kernel(u, True), 0, hi, epsabs=1e-13, epsrel=1e-12)[0]
    log_z = np.log(2) + d / 2 * np.log(np.pi) - gammaln(d / 2) - d / 2 * np.log(mu) + np.log(z)
    entropy = np.log(2) + (d + 1) / 2 * np.log(np.pi) - gammaln((d + 1) / 2)
    distortion = m / z / mu
    rate = entropy - log_z - mu * distortion
    residual = d / 2 * np.log(2 * np.pi * np.e * distortion / d) - log_z - mu * distortion
    return {"dimension": d, "criterion": criterion, "mu": mu,
            "D": float(distortion), "R_nats": float(rate),
            "residual_nats": float(residual), "residual_over_D": float(residual / distortion)}


def surface_geometry(kind, n, major=2.0, minor=1.0, a=1.3, c=0.7, kappa=0.0):
    """Positive meridian quadrature weights, including area and density."""
    if kind == "torus":
        t = 2 * np.pi * np.arange(n) / n
        wt = np.full(n, 2 * np.pi / n)
        f, g = major + minor * np.cos(t), minor * np.sin(t)
        fp, gp = -minor * np.sin(t), minor * np.cos(t)
        fpp, gpp = -minor * np.cos(t), -minor * np.sin(t)
    elif kind == "ellipsoid":
        u, wt = roots_legendre(n)
        t, wt = (u + 1) * np.pi / 2, wt * np.pi / 2
        f, g = a * np.sin(t), c * np.cos(t)
        fp, gp = a * np.cos(t), -c * np.sin(t)
        fpp, gpp = -a * np.sin(t), -c * np.cos(t)
    else:
        raise ValueError(kind)
    speed = np.hypot(fp, gp)
    area_w = 2 * np.pi * f * speed * wt
    if np.any(area_w <= 0):
        raise ValueError("surface quadrature must have positive area weights")
    density = np.exp(kappa * np.cos(t))
    normalizer = float(area_w @ density)
    p = density / normalizer
    mass = p * area_w
    k1 = (fp * gpp - gp * fpp) / speed ** 3
    k2 = gp / (f * speed)
    return {"t": t, "f": f, "g": g, "mass": mass,
            "entropy": float(-mass @ np.log(p)),
            "coefficient": float(-mass @ (k1 - k2) ** 2 / 16),
            "local_A": (k1 - k2) ** 2 / 16}


def surface_ba(geom, mu, tolerance=1e-12, max_iter=50000):
    """Symmetry-reduced BA with the reproduction alphabet on the surface.

    K_ij=exp(-mu*e_ij) I0e(2 mu f_i f_j), where e_ij is the
    squared meridian distance. I0e already includes exp(-2 mu f_i f_j).
    At any iterate, log(max_j c_j) bounds the Lagrangian optimality gap.
    The reported mutual information uses the induced output marginal.
    """
    f, g, w = geom["f"], geom["g"], geom["mass"]
    e = (f[:, None] - f[None, :]) ** 2 + (g[:, None] - g[None, :]) ** 2
    aa = 2 * f[:, None] * f[None, :]
    ex = np.exp(-mu * e)
    bessel0, bessel1 = i0e(mu * aa), i1e(mu * aa)
    kernel = ex * bessel0
    moment = kernel * e + ex * aa * (bessel0 - bessel1)
    q = w.copy()
    for iteration in range(max_iter + 1):
        zx = kernel @ q
        if np.any(zx <= 0) or not np.all(np.isfinite(zx)):
            raise FloatingPointError("unresolved Gibbs kernel")
        cvec = kernel.T @ (w / zx)
        gap = max(0.0, float(np.log(cvec.max())))
        if gap <= tolerance:
            break
        q *= cvec
        q /= q.sum()
    else:
        raise RuntimeError(f"BA did not converge: mu={mu}, gap={gap:g}")
    output_mass = q * cvec
    output_mass /= output_mass.sum()
    kl = float(np.sum(output_mass * np.log(np.maximum(cvec, np.finfo(float).tiny))))
    distortion = float(w @ ((moment @ q) / zx))
    cross_lagrangian = float(-w @ np.log(zx))
    rate = cross_lagrangian - mu * distortion - kl
    residual = rate - (geom["entropy"] - np.log(np.pi * np.e * distortion))
    local_d = ((w / zx) @ moment) / cvec
    return {"mu": mu, "D": distortion, "R_nats": rate,
            "residual_nats": residual, "residual_over_D": residual / distortion,
            "lagrangian_gap_nats": gap, "iterations": iteration,
            "output_marginal_error_l1": float(np.abs(output_mass - q).sum()),
            "local_scaled_D": (mu * mu * (local_d - 1 / mu)).tolist()}


def projection_sphere(d, variance, order=80):
    """Exact Gaussian expectation of geodesic squared error on S^d(1).

    The radial normal coordinate is N(0,v); tangent squared norm/(2v)
    has Gamma(d/2,1) law. Generalized Laguerre quadrature is essential
    for half-integer d/2, in particular for the circle.
    """
    tangential, wt = roots_genlaguerre(order, d / 2 - 1)
    z, wz = roots_hermitenorm(order)
    wt = wt / np.exp(gammaln(d / 2))
    wz = wz / np.sqrt(2 * np.pi)
    angle = np.arctan2(np.sqrt(2 * variance * tangential[None, :]),
                      1 + np.sqrt(variance) * z[:, None])
    return float(wz @ (angle ** 2) @ wt)


def outage_two_circles(sigma_axis, samples, seed):
    """Two equiprobable circles; sigma_axis is per-coordinate noise std."""
    r0, r1, da = 1.0, 1.6, 32
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, 2, size=samples)
    radii = np.where(labels == 0, r0, r1)
    # Rotational symmetry fixes each noiseless point to (r,0).
    y = np.column_stack([radii, np.zeros(samples)]) + rng.normal(scale=sigma_axis, size=(samples, 2))
    predictions = (np.linalg.norm(y, axis=1) > (r0 + r1) / 2).astype(int)
    failures = int(np.count_nonzero(labels != predictions))
    phat = failures / samples
    # Wilson binomial interval includes zero-event and rare-event runs.
    z = norm.ppf(0.975)
    den = 1 + z * z / samples
    center = (phat + z * z / (2 * samples)) / den
    half = z / den * np.sqrt(phat * (1 - phat) / samples + z * z / (4 * samples * samples))
    threshold = ((r0 + r1) / 2 / sigma_axis) ** 2
    exact = (ncx2.sf(threshold, 2, (r0 / sigma_axis) ** 2)
             + ncx2.cdf(threshold, 2, (r1 / sigma_axis) ** 2)) / 2
    return {"sigma_axis": sigma_axis, "sigma_total": sigma_axis * np.sqrt(da),
            "ambient_dimension": da, "samples": samples, "seed": seed,
            "failures": failures, "p_out_mc": phat, "ci95_low": center - half,
            "ci95_high": center + half, "p_out_exact": float(exact),
            "normal_tail_approx": float(norm.sf((r1 - r0) / (2 * sigma_axis)))}


def write_csv(path, rows):
    fields = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plots(result, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.labelsize": 8, "legend.fontsize": 7,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.8, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "legend.frameon": False})
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#555555"]
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 2.7), layout="constrained")
    for criterion, color, marker, label in [("geodesic", colors[0], "o", "Geodesic"),
                                           ("chordal", colors[1], "s", "Chordal")]:
        rows = [r for r in result["sphere_coefficients"] if r["criterion"] == criterion]
        axes[0].plot([r["dimension"] for r in rows], [r["predicted"] for r in rows],
                     color=color, linestyle="-" if criterion == "geodesic" else "--", label=label + " theory")
        axes[0].plot([r["dimension"] for r in rows], [r["fitted"] for r in rows],
                     color=color, marker=marker, fillstyle="none", linestyle="none", label=label + " integral")
    axes[0].set(xlabel="Intrinsic dimension d", ylabel="RD coefficient (nats / squared length)", title="(a) Uniform spheres")
    axes[0].axhline(0, color="0.7", lw=.6); axes[0].legend(loc="upper left")
    for i, surface in enumerate(result["surfaces"]):
        rows = surface["curve"]
        axes[1].plot([r["D"] for r in rows], [r["residual_over_D"] for r in rows],
                     marker=["o", "s", "^", "D", "v"][i], color=colors[i],
                     label=surface["label"], markersize=4)
        axes[1].scatter([0], [surface["predicted"]], marker="x", color=colors[i], s=30)
    axes[1].set(xlabel="Distortion D (squared length)", ylabel=r"$[R(D)-R_0(D)]/D$ (nats / squared length)", title="(b) Chordal distortion on surfaces")
    axes[1].legend(loc="center right", bbox_to_anchor=(1, .40)); axes[1].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    for ax in axes: ax.grid(alpha=.18)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / f"rd_validation.{ext}", dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 2.7), layout="constrained")
    for i, d in enumerate((1, 2, 3, 4, 5)):
        rows = [r for r in result["projection_curve"] if r["dimension"] == d]
        axes[0].plot([r["variance_axis"] for r in rows], [r["second_order_ratio"] for r in rows],
                     color=colors[i], marker=["o", "s", "^", "D", "v"][i], markersize=4, label=f"d = {d}")
        axes[0].scatter([0], [d * (5 - 2 * d) / 3], marker="x", color=colors[i], s=30)
    axes[0].set(xlabel="Per-axis noise variance v", ylabel=r"$[\mathbb{E}d_g^2-dv]/v^2$", title="(a) Projection on unit spheres")
    axes[0].legend(ncol=2, loc="center right", bbox_to_anchor=(1, .35)); axes[0].axhline(0, color="0.7", lw=.6)
    rows = result["outage"]
    x = [r["sigma_axis"] for r in rows]
    y = np.array([r["p_out_mc"] for r in rows])
    axes[1].errorbar(x, y, yerr=np.array([y - [r["ci95_low"] for r in rows],
                                        [r["ci95_high"] for r in rows] - y]),
                     fmt="o", color=colors[0], capsize=2, markersize=4, label="Monte Carlo, 95% CI")
    axes[1].plot(x, [r["p_out_exact"] for r in rows], color=colors[0], label="Exact radial law")
    axes[1].plot(x, [r["normal_tail_approx"] for r in rows], "--", color=colors[1], label=r"$Q(\Delta/(2\sqrt{v}))$")
    axes[1].set(xlabel=r"Per-axis standard deviation $\sqrt{v}$", ylabel="Semantic outage probability", yscale="log", title="(b) Two circles, separation 0.6")
    axes[1].legend(loc="lower right")
    for ax in axes: ax.grid(alpha=.18)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / f"channel_validation.{ext}", dpi=220)
    plt.close(fig)


def run(output, n_meridian=768, samples=400000, seed=20261019):
    started = time.time()
    output.mkdir(parents=True, exist_ok=True)
    result = {"metadata": {"schema_version": 1, "seed": seed,
              "python": sys.version.split()[0], "numpy": np.__version__,
              "scipy": scipy.__version__, "platform": platform.platform(),
              "rate_unit": "nat", "length_unit": "sphere radius or minor torus radius",
              "n_meridian": n_meridian, "ba_tolerance_nats": {"torus": 1e-12, "ellipsoid": 1e-9},
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "command": " ".join(sys.argv)},
              "sphere_curve": [], "sphere_coefficients": [], "surfaces": [],
              "projection_curve": [], "projection_coefficients": [], "outage": []}
    for criterion in ("geodesic", "chordal"):
        for d in (1, 2, 3, 5, 8):
            rows = [sphere_rd(d, mu, criterion) for mu in (500., 1000., 2000., 4000.)]
            fitted = extrapolate([r["D"] for r in rows], [r["residual_over_D"] for r in rows])
            prediction = (d - 1) / 6 if criterion == "geodesic" else (d - 2) / 8
            result["sphere_curve"].extend(rows)
            result["sphere_coefficients"].append({"dimension": d, "criterion": criterion,
                "predicted": prediction, "fitted": fitted, "absolute_error": abs(fitted - prediction)})
    print("Sphere RD integrals completed.", flush=True)
    specs = [("Torus R=2, k=0", {"kind": "torus", "major": 2., "kappa": 0.}),
             ("Torus R=2, k=1", {"kind": "torus", "major": 2., "kappa": 1.}),
             ("Torus R=2, k=-1", {"kind": "torus", "major": 2., "kappa": -1.}),
             ("Torus R=3, k=0", {"kind": "torus", "major": 3., "kappa": 0.}),
             ("Ellipsoid 1.3, 0.7", {"kind": "ellipsoid", "a": 1.3, "c": .7})]
    for label, spec in specs:
        geom = surface_geometry(n=n_meridian, **spec)
        mus = (100., 200., 400., 800.) if spec["kind"] == "torus" else (200., 400., 800., 1600.)
        tol = 1e-12 if spec["kind"] == "torus" else 1e-9
        curve = [surface_ba(geom, mu, tolerance=tol) for mu in mus]
        fitted = extrapolate([r["D"] for r in curve], [r["residual_over_D"] for r in curve])
        surface = {"label": label, "parameters": spec, "predicted": geom["coefficient"],
                   "fitted": fitted, "absolute_error": abs(fitted - geom["coefficient"]), "curve": curve}
        if label == "Torus R=2, k=0":
            refined = surface_geometry(n=2*n_meridian, **spec)
            fine = surface_ba(refined, 800.)
            surface["grid_refinement"] = {"n_coarse": n_meridian, "n_fine": 2*n_meridian,
                "mu": 800., "rate_difference_nats": abs(fine["R_nats"] - curve[-1]["R_nats"]),
                "coefficient_ratio_difference": abs(fine["residual_over_D"] - curve[-1]["residual_over_D"])}
        # The full profile is reproducible; retaining only selected locations
        # avoids a large array in every aggregate result row.
        if label == "Torus R=2, k=0":
            ii = [0, n_meridian//4, n_meridian//2]
            surface["local_profile"] = {"theta": geom["t"][ii].tolist(),
                "measured": np.asarray(curve[-1]["local_scaled_D"])[ii].tolist(),
                "leading_prediction": geom["local_A"][ii].tolist(),
                "conditioning": "reverse channel: source distortion conditional on reproduction ring"}
        for row in curve: row.pop("local_scaled_D")
        result["surfaces"].append(surface)
        print(f"{label}: predicted={geom['coefficient']:+.8f}, fitted={fitted:+.8f}", flush=True)
    for d in (1, 2, 3, 4, 5):
        rows = []
        for v in (.001, .002, .004, .008):
            mse = projection_sphere(d, v)
            rows.append({"dimension": d, "variance_axis": v, "geodesic_mse": mse,
                "second_order_ratio": (mse - d*v)/(v*v)})
        fitted = extrapolate([r["variance_axis"] for r in rows], [r["second_order_ratio"] for r in rows])
        predicted = d*(5-2*d)/3
        result["projection_curve"].extend(rows)
        result["projection_coefficients"].append({"dimension": d, "predicted": predicted,
            "fitted": fitted, "absolute_error": abs(fitted-predicted),
            "quadrature_refinement_abs_mse": abs(projection_sphere(d,.008,120)-rows[-1]["geodesic_mse"])})
    for i, s in enumerate((.08, .10, .12, .15, .177, .20, .25)):
        result["outage"].append(outage_two_circles(s, samples, seed+i))
    result["checks"] = {
        "sphere_coefficient_max_absolute_error": max(r["absolute_error"] for r in result["sphere_coefficients"]),
        "surface_coefficient_max_absolute_error": max(r["absolute_error"] for r in result["surfaces"]),
        "projection_coefficient_max_absolute_error": max(r["absolute_error"] for r in result["projection_coefficients"]),
        "maximum_ba_gap_nats": max(r["lagrangian_gap_nats"] for s in result["surfaces"] for r in s["curve"])}
    result["metadata"]["elapsed_seconds"] = time.time()-started
    (output/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    for name in ("sphere_curve", "sphere_coefficients", "projection_curve", "projection_coefficients", "outage"):
        write_csv(output/(name+".csv"), result[name])
    write_csv(output/"surface_coefficients.csv", [{k:s[k] for k in ("label","predicted","fitted","absolute_error")} for s in result["surfaces"]])
    plots(result, output)
    print(json.dumps(result["checks"], indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1]/"results"/"paper1")
    parser.add_argument("--n-meridian", type=int, default=768)
    parser.add_argument("--samples", type=int, default=400000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        run(args.output, args.n_meridian, args.samples, args.seed)
