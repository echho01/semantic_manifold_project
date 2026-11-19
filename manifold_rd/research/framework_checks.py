"""Numerical checks of the manifold framework: generalized, knowledge-conditioned
and task-oriented semantic rate-distortion, and channel-induced semantic distortion.

  homogeneous  exact parametric RDFs of S^1, S^2, S^3 and the flat torus T^2 under squared
               geodesic distortion: slope m/2 and curvature coefficient Scal/(6m);
  nonuniform   a tilted zonal density on S^2: achievable kernel and dual bound against the
               same coefficient, which contains no Fisher-information term;
  spheroid     a surface of varying curvature with a nonuniform density: the coefficient is the
               density-weighted mean curvature, not the area average;
  knowledge    a circle whose shared knowledge K changes the source law, the semantic
               metric, or both; the gain R_0(D) - R_K(D) against I(U;K) - E ln w_K;
  task         a torus (theta, phi) whose task posterior depends on theta only; generalized,
               task-oriented and joint (KL + eta * distance^2) RDFs and their asymptotes;
  channel      a one-shot analog circle map over two real Gaussian uses: exact semantic
               distortion against the Fisher-information prediction.

Run: python framework_checks.py [--output results] [--figure ../figures/framework.png]
Finite-grid values are numerical evaluations of achievable Blahut points, not certificates.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.special import ndtr

from manifold_rd import finite_ba, posterior_kl_bits
from theory_checks import LN2, circle, fisher_speed_and_curvature

np.seterr(divide="ignore", over="ignore", invalid="ignore")


# --------------------------------------------------------------------------
# A. Homogeneous manifolds: exact RDF, slope and curvature coefficient
# --------------------------------------------------------------------------
SPACES = {
    # name: (dimension, volume, scalar curvature, radial density of the volume element on [0, pi])
    "S1": (1, 2 * np.pi, 0.0, lambda x: 2.0 + 0 * x),
    "S2": (2, 4 * np.pi, 2.0, lambda x: 2 * np.pi * np.sin(x)),
    "S3": (3, 2 * np.pi**2, 6.0, lambda x: 4 * np.pi * np.sin(x) ** 2),
}


def radial_moments(density, s):
    """Z_s and the mean squared distance of the Gibbs law exp(-s x^2) on a radial model."""
    cut = min(np.pi, 12.0 / np.sqrt(s))
    z = integrate.quad(lambda x: np.exp(-s * x * x) * density(x), 0, cut, epsabs=0, epsrel=1e-13)[0]
    d = integrate.quad(lambda x: x * x * np.exp(-s * x * x) * density(x), 0, cut, epsabs=0, epsrel=1e-13)[0]
    return z, d / z


def homogeneous_point(name, s):
    """Exact (distortion, rate in nats) at Gibbs scale s; T2 is a product of two circles."""
    if name == "T2":
        z, d = radial_moments(SPACES["S1"][3], s)
        return 2 * d, 2 * (np.log(2 * np.pi) - np.log(z) - s * d)
    m, vol, _, density = SPACES[name]
    z, d = radial_moments(density, s)
    return d, np.log(vol) - np.log(z) - s * d


def homogeneous():
    info = dict(SPACES, T2=(2, 4 * np.pi**2, 0.0, None))
    rows, curves = [], {}
    for name, (m, vol, scal, _) in info.items():
        curve = []
        for s in np.geomspace(0.5, 2e4, 60):
            d, rate = homogeneous_point(name, s)
            lead = np.log(vol) - 0.5 * m * np.log(2 * np.pi * np.e * d / m)
            curve.append((d, rate / LN2, (rate - lead) / d))
        curves[name] = curve
        for target in (0.1, 0.01, 0.001):
            near = min(curve, key=lambda c: abs(np.log(c[0] / target)))
            rows.append(dict(space=name, dimension=m, scalar_curvature=scal, distortion=near[0],
                             rate_bits=near[1], coefficient=near[2], predicted=scal / (6 * m)))
        # Local slope dR/dln(1/D) in nats at small distortion.
        (d1, r1, _), (d2, r2, _) = curve[-1], curve[-6]
        rows[-1]["slope_nats"] = (r1 - r2) * LN2 / np.log(d2 / d1)
    return rows, curves


def nonuniform_sphere(tilt, t, nodes=4000):
    """Zonal source p = (1 + tilt cos(theta)) / 4pi on S^2 under squared geodesic distortion.

    Returns the exact distortion and rate (nats) of the source-weighted Gibbs kernel
    p(y) exp(-t d^2) / (T_t p)(x), the dual lower bound h(p) - ln(T_t 1) - t D at the same
    distortion, the entropy h(p) and the Fisher information E|grad ln p|^2. Zonal functions are
    convolved with the kernel through its Legendre coefficients (Funk-Hecke).
    """
    from numpy.polynomial import legendre
    lmax = int(12 * np.sqrt(t)) + 20
    x, w = legendre.leggauss(nodes)                      # x = cos(theta)
    gamma = np.linspace(0, np.pi, 200001)
    step = np.full(len(gamma), gamma[1])
    step[[0, -1]] /= 2
    basis = legendre.legvander(np.cos(gamma), lmax)
    weight = np.exp(-t * gamma**2) * np.sin(gamma) * step
    lam, mu = 2 * np.pi * (weight @ basis), 2 * np.pi * ((gamma**2 * weight) @ basis)
    nodes_basis = legendre.legvander(x, lmax)

    def convolve(coefficients, f):
        return nodes_basis @ (coefficients * (np.arange(lmax + 1) + 0.5) * ((f * w) @ nodes_basis))

    mean = lambda f: 2 * np.pi * np.sum(w * f)           # integral of a zonal function over S^2
    p = (1 + tilt * x) / (4 * np.pi)
    smooth = convolve(lam, p)
    distortion = mean(p * convolve(mu, p) / smooth)
    output_ratio = convolve(lam, p / smooth)             # output density divided by p
    rate = -t * distortion - mean(p * np.log(smooth)) - mean(p * output_ratio * np.log(output_ratio))
    entropy = -mean(p * np.log(p))
    return dict(tilt=tilt, scale=t, distortion=distortion, rate_nats=rate,
                dual_lower_nats=entropy - np.log(lam[0]) - t * distortion, entropy_nats=entropy,
                fisher_information=mean(p * tilt**2 * (1 - x * x) / (1 + tilt * x) ** 2))


def spheroid(a, c, tilt, t, nz=64, nphi=64, nr=600, cut=7.0):
    """Spheroid x^2/a^2 + y^2/a^2 + z^2/c^2 = 1 with density p ~ 1 + tilt (z/c)^2, squared geodesic cost.

    The Gaussian curvature varies with latitude, so the density-weighted mean E_p[Scal] differs
    from the area average. Around each base point the kernel is integrated in geodesic polar
    coordinates: geodesics are shot by RK4 in the embedding and the volume element is the Jacobi
    field J'' = -K J. Returns normalized remainders of the source-weighted Gibbs kernel (upper)
    and of a dual bound (lower), to be compared with E_p[Scal]/(6m) = E_p[K]/6.
    """
    from numpy.polynomial import legendre as L
    x, w = L.leggauss(nz); z = c * x
    rho = a * np.sqrt(1 - x * x)
    area = 2 * np.pi * np.sqrt(rho**2 + (a * a * z / c**2) ** 2) * c * w      # zonal area weights
    K = lambda P: 1.0 / (a**4 * c**2 * ((P[..., 0]**2 + P[..., 1]**2) / a**4 + P[..., 2]**2 / c**4) ** 2)
    grad = lambda P: P * np.array([2 / a**2, 2 / a**2, 2 / c**2])
    hdiag = np.array([2 / a**2, 2 / a**2, 2 / c**2])
    def rhs(P, V, J, dJ):
        g = grad(P)
        acc = -((V * V * hdiag).sum(-1) / (g * g).sum(-1))[..., None] * g
        return V, acc, dJ, -K(P) * J
    # base points and tangent frames
    P0 = np.stack([rho, 0 * rho, z], -1)
    e1 = np.stack([-a * a * z / c**2, 0 * z, rho], -1); e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)   # meridian tangent (up to sign)
    e2 = np.array([0.0, 1.0, 0.0])
    phi = 2 * np.pi * np.arange(nphi) / nphi
    V = np.cos(phi)[None, :, None] * e1[:, None, :] + np.sin(phi)[None, :, None] * e2
    P = np.repeat(P0[:, None, :], nphi, 1); J = np.zeros((nz, nphi)); dJ = np.ones((nz, nphi))
    h = cut / np.sqrt(t) / nr
    Z = np.empty((nr + 1, nz, nphi)); Js = np.empty_like(Z); Z[0] = P[..., 2]; Js[0] = 0
    for i in range(nr):
        k1 = rhs(P, V, J, dJ)
        k2 = rhs(P + h/2*k1[0], V + h/2*k1[1], J + h/2*k1[2], dJ + h/2*k1[3])
        k3 = rhs(P + h/2*k2[0], V + h/2*k2[1], J + h/2*k2[2], dJ + h/2*k2[3])
        k4 = rhs(P + h*k3[0], V + h*k3[1], J + h*k3[2], dJ + h*k3[3])
        P = P + h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); V = V + h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        J = J + h/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]); dJ = dJ + h/6*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
        Z[i+1] = P[..., 2]; Js[i+1] = J
    assert Js.min() >= 0, "conjugate point inside the cut radius"
    r = h * np.arange(nr + 1)
    simpson = np.ones(nr + 1); simpson[1:-1:2] = 4; simpson[2:-1:2] = 2; simpson *= h / 3
    base = (simpson * np.exp(-t * r * r))[:, None, None] * Js * (2 * np.pi / nphi)   # kernel * volume element
    Vand = L.legvander(x, nz - 1)
    def at(f_nodes):                       # zonal function given at nodes, evaluated at all shot points
        coef = np.linalg.solve(Vand, f_nodes)
        return L.legval(np.clip(Z / c, -1, 1), coef)
    T = lambda f, extra=1.0: (base * extra * at(f)).sum((0, 2))
    p = 1 + tilt * x * x; p /= (p * area).sum()
    Kn = K(P0); E = lambda f: float((p * area * f).sum())
    g = T(p)
    D = E(T(p, (r * r)[:, None, None]) / g)
    ratio = T(p / g)
    I = -t * D - E(np.log(g)) - E(ratio * np.log(ratio))
    hent = -E(np.log(p))
    wdual = t / np.pi * T(np.exp(2 * Kn / (12 * t)))
    low = hent + np.log(t / np.pi) + E(2 * Kn) / (12 * t) - np.log(wdual.max()) - t * D
    lead = hent - np.log(np.pi * np.e * D)
    return dict(a=a, c=c, tilt=tilt, scale=t, distortion=D, upper_coefficient=float((I - lead) / D),
                lower_coefficient=float((low - lead) / D), predicted=E(Kn) / 6,
                area_average=float((Kn * area).sum() / area.sum()) / 6,
                gauss_bonnet=float((Kn * area).sum() / (4 * np.pi)))


def variable_curvature():
    return [spheroid(1.0, 0.7, tilt, t) for tilt in (3.0, -0.9) for t in (30, 100, 300, 1000)]


def nonuniform():
    """Curvature coefficient for nonuniform densities: both bounds against Scal/(6m) = 1/6."""
    rows = []
    for tilt in (0.0, 0.5, 0.8):
        for t in (10, 30, 100, 300):
            r = nonuniform_sphere(tilt, t)
            lead = r["entropy_nats"] - np.log(np.pi * np.e * r["distortion"])
            r.update(upper_coefficient=(r["rate_nats"] - lead) / r["distortion"],
                     lower_coefficient=(r["dual_lower_nats"] - lead) / r["distortion"],
                     gap_over_d2=(r["rate_nats"] - r["dual_lower_nats"]) / r["distortion"] ** 2,
                     predicted=1 / 6,
                     # What the coefficient would be if a Fisher term J/(4t) = J D/4 survived.
                     coefficient_with_fisher_term=1 / 6 + r["fisher_information"] / 4)
            rows.append(r)
    return rows


# --------------------------------------------------------------------------
# B. Knowledge: law, metric, or both
# --------------------------------------------------------------------------
def rd_curve(mass, cost, slopes):
    pts = [finite_ba(mass, cost, b, tol=1e-9, max_iter=20000) for b in slopes]
    return np.array([(p.distortion, p.rate) for p in pts])


def rate_at(curve, target):
    order = np.argsort(curve[:, 0])
    return float(np.interp(np.log(target), np.log(curve[order, 0]), curve[order, 1]))


def knowledge(n=720):
    theta = 2 * np.pi * np.arange(n) / n
    slopes = np.geomspace(0.5, 1500, 34)

    def cost(scale, tilt):               # squared geodesic distance for the metric (scale (1 + tilt cos))^2 dtheta^2
        arc = scale * (theta + tilt * np.sin(theta))
        gap = np.abs(arc[:, None] - arc)
        return np.minimum(gap, 2 * np.pi * scale - gap) ** 2

    uniform = np.full(n, 1 / n)
    tilted = (1 + 0.8 * np.cos(theta)) / n
    base = rd_curve(uniform, cost(1.0, 0.0), slopes)
    mutual = float(np.log2(n) + (tilted * np.log2(tilted)).sum())           # I(U;K), bits, mirror-symmetric K
    cases = [
        ("law only", tilted, 1.0, 0.0),
        ("metric only", uniform, 1.0, 0.5),
        ("coarser where likely", tilted, 1.0, -0.5),
        ("finer where likely", tilted, 1.25, 0.5),
    ]
    rows = []
    for name, mass, scale, tilt in cases:
        curve = rd_curve(mass, cost(scale, tilt), slopes)
        log_w = np.log2(scale * (1 + tilt * np.cos(theta)))
        info = mutual if mass is tilted else 0.0
        predicted = info - float(mass @ log_w)
        row = dict(case=name, mutual_information_bits=info, mean_log_metric_bits=float(mass @ log_w),
                   predicted_gain_bits=predicted)
        for target in (0.3, 0.1, 0.03, 0.01, 0.003):
            row[f"gain_at_{target}"] = rate_at(base, target) - rate_at(curve, target)
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# C. Task: a torus whose posterior depends on one coordinate
# --------------------------------------------------------------------------
def circle_gibbs(s):
    """Exact (distortion, rate in nats) of the uniform circle under cost x^2 at Gibbs scale s."""
    z, d = radial_moments(SPACES["S1"][3], s)
    return d, np.log(2 * np.pi) - np.log(z) - s * d


def task(n=720, amplitude=0.2, eta=0.01):
    theta = 2 * np.pi * np.arange(n) / n
    q = circle(theta, amplitude)
    p = np.full(n, 1 / n)
    kl = posterior_kl_bits(q, q) * LN2                              # nats
    gap = np.abs(theta[:, None] - theta)
    arc2 = np.minimum(gap, 2 * np.pi - gap) ** 2
    off = ~np.eye(n, dtype=bool)
    lipschitz = float((kl[off] / arc2[off]).max())                  # KL <= L^2 arc^2 on the grid
    g, _ = fisher_speed_and_curvature(theta, amplitude)
    h_fisher = float(np.log(2 * np.pi) + 0.5 * np.mean(np.log(g)))
    h_joint = float(np.log(4 * np.pi**2) + 0.5 * np.mean(np.log(0.5 * g + eta)) + 0.5 * np.log(eta))
    slopes = np.geomspace(3, 6000, 40)                              # bits per nat of distortion
    task_curve, joint_curve, general_curve = [], [], []
    for beta in slopes:
        a = finite_ba(p, kl, beta, tol=1e-9, max_iter=20000)
        task_curve.append((a.distortion, a.rate))
        b = finite_ba(p, kl + eta * arc2, beta, tol=1e-9, max_iter=20000)
        d_phi, r_phi = circle_gibbs(beta * LN2 * eta)               # free coordinate, cost eta * arc^2
        joint_curve.append((b.distortion + eta * d_phi, b.rate + r_phi / LN2))
        scale = lipschitz + eta                                     # calibrated generalized cost on the torus
        d0, r0 = circle_gibbs(beta * LN2 * scale)
        general_curve.append((2 * scale * d0, 2 * r0 / LN2))
    task_curve, joint_curve, general_curve = map(np.array, (task_curve, joint_curve, general_curve))
    rows = []
    for target in (0.1, 0.03, 0.01, 0.003, 0.001, 0.0003):
        r_task, r_joint, r_gen = (rate_at(c, target) for c in (task_curve, joint_curve, general_curve))
        rows.append(dict(distortion_nats=target, task_bits=r_task, joint_bits=r_joint, generalized_bits=r_gen,
                         task_asymptote_bits=(h_fisher - 0.5 * np.log(4 * np.pi * np.e * target)) / LN2,
                         joint_asymptote_bits=(h_joint - np.log(np.pi * np.e * target)) / LN2))
    summary = dict(grid=n, eta=eta, lipschitz_sq=lipschitz, h_fisher_nat=h_fisher, h_joint_nat=h_joint,
                   half_mean_fisher_speed_sq=0.5 * float(np.mean(g)))
    curves = dict(task=task_curve.tolist(), joint=joint_curve.tolist(), generalized=general_curve.tolist())
    return summary, rows, curves


# --------------------------------------------------------------------------
# D. Channel-induced semantic distortion of a one-shot analog map
# --------------------------------------------------------------------------
def channel(amplitude=0.2):
    """Send rho (cos theta, sin theta) over two unit-variance real uses; decode the angle."""
    theta = 2 * np.pi * (np.arange(360) + 0.5) / 360
    g, _ = fisher_speed_and_curvature(theta, amplitude)
    q = circle(theta, amplitude)
    rows = []
    for rho in (1.0, 2.0, 4.0, 8.0, 16.0):
        def density(phi):
            a = rho * np.cos(phi)
            return (np.exp(-rho**2 / 2) / (2 * np.pi)
                    + a * np.exp((a * a - rho * rho) / 2) * ndtr(a) / np.sqrt(2 * np.pi))
        nodes, weights = np.polynomial.legendre.leggauss(4000)
        phi, wts = np.pi * nodes, np.pi * weights                   # the whole circle of angular errors
        f = density(phi)
        if abs(float(f @ wts) - 1) > 1e-10:
            raise ArithmeticError("angular error density does not integrate to one")
        second = float((phi**2 * f) @ wts)
        v = 1 / rho**2                                              # tr(G G_obs^{-1}) with G = 1
        kl = np.array([(q * (np.log(q) - np.log(circle(theta + e, amplitude)))).sum(1).mean() for e in phi])
        rows.append(dict(rho=rho, fisher_bound=v, second_order=v + v * v, exact_sq_error=second,
                         ratio_to_fisher=second / v, task_kl_nats=float(kl @ (f * wts)),
                         task_prediction_nats=0.5 * float(np.mean(g)) * v))
    return rows


def plot(hom_curves, task_summary, task_curves, path):
    import matplotlib.pyplot as plt

    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(7.2, 2.5), layout="constrained")
    colors = {"S1": "#6a8f3c", "S2": "#124971", "S3": "#a32437", "T2": "#c98a1c"}
    labels = {"S1": r"$S^1$, $m=1$", "S2": r"$S^2$, $m=2$", "S3": r"$S^3$, $m=3$", "T2": r"$T^2$, $m=2$"}
    for name, curve in hom_curves.items():
        c = np.array(curve)
        ax.semilogx(c[:, 0], c[:, 1], color=colors[name], label=labels[name], lw=1.2)
        bx.semilogx(c[:, 0], c[:, 2], color=colors[name], lw=1.2)
    ax.set(xlabel="Distortion $D$ (squared distance)", ylabel="Rate (bits)", xlim=(1e-3, 1), ylim=(0, 16),
           title="(a) Dimension sets the slope")
    ax.legend(frameon=False, fontsize=6.5)
    for level in (1 / 6, 1 / 3):
        bx.axhline(level, color=".5", ls=":", lw=.8)
    bx.set(xlabel="Distortion $D$", ylabel=r"$(R-R_{\rm lead})/D$ (nats)", xlim=(1e-3, 1), ylim=(-0.05, 0.45),
           title=r"(b) Curvature term ${\rm Scal}/(6m)$")
    t, j, g0 = (np.array(task_curves[k]) for k in ("task", "joint", "generalized"))
    d = np.geomspace(2e-4, 0.3, 100)
    cx.semilogx(g0[:, 0], g0[:, 1], color="#c98a1c", lw=1.2, label="Generalized (calibrated)")
    cx.semilogx(j[:, 0], j[:, 1], color="#124971", lw=1.2, label=r"Joint, $\eta=0.01$")
    cx.semilogx(t[:, 0], t[:, 1], color="#a32437", lw=1.2, label="Task")
    cx.semilogx(d, (task_summary["h_fisher_nat"] - 0.5 * np.log(4 * np.pi * np.e * d)) / LN2, ":", color="#a32437", lw=.8)
    cx.semilogx(d, (task_summary["h_joint_nat"] - np.log(np.pi * np.e * d)) / LN2, ":", color="#124971", lw=.8)
    cx.set(xlabel="Distortion $D$ (nats)", ylabel="Rate (bits)", xlim=(2e-4, 0.3), ylim=(0, 12),
           title="(c) Task removes dimensions")
    cx.legend(frameon=False, fontsize=6.5)
    for axis in (ax, bx, cx):
        axis.grid(alpha=.2, which="both")
        axis.title.set_fontsize(8.5)
        axis.tick_params(labelsize=7)
        axis.xaxis.label.set_fontsize(8)
        axis.yaxis.label.set_fontsize(8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--knowledge-grid", type=int, default=720)
    args = parser.parse_args()
    hom_rows, hom_curves = homogeneous()
    task_summary, task_rows, task_curves = task()
    record = dict(homogeneous=hom_rows, nonuniform_sphere=nonuniform(), variable_curvature=variable_curvature(), knowledge=knowledge(args.knowledge_grid), task_summary=task_summary,
                  task=task_rows, task_curves=task_curves, channel=channel())
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "framework_checks.json").write_text(json.dumps(record, indent=1))
    if args.figure:
        plot(hom_curves, task_summary, task_curves, args.figure)
    print(f"wrote {args.output / 'framework_checks.json'}")


if __name__ == "__main__":
    main()
