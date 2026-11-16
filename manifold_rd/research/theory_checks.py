"""Numerical checks for the structural results added to the Transactions manuscript.

Everything here is a deterministic computation on the three-class posterior circle
q_i(theta) = 1/3 + a cos(theta - 2 pi i / 3).  Four blocks are produced:

  budget   exact loss budget of a complete link (representation, compression,
           channel, decoder) against the directly evaluated log loss;
  gap      object versus predictive receivers at equal information on dense
           alphabets, with the curvature and reach bounds of the manuscript;
  highres  small-distortion behaviour against the task-rank asymptote and the
           mixture-curvature constant;
  opta     distortion floors implied by the certified finite bounds at the
           capacities of the two-use Gaussian link, plus the knowledge saving.

Run: python theory_checks.py [--output results] [--figure ../figures/gap_highres.png]
Dense-alphabet values are numerical evaluations, not certificates.  The 24-state
quantities reuse the certified solver in manifold_rd.py.
"""
import argparse
import json
from math import erfc
from pathlib import Path

import numpy as np

from manifold_rd import conditional_ba, finite_ba, posterior_kl_bits

LN2 = np.log(2.0)
# Accelerate-backed matmul raises spurious floating-point flags on some macOS builds.
np.seterr(divide="ignore", over="ignore", invalid="ignore")
PHASES = 2 * np.pi * np.arange(3) / 3
CENTER = np.full(3, 1 / 3)


def circle(angle, amplitude=0.2):
    return 1 / 3 + amplitude * np.cos(np.asarray(angle)[..., None] - PHASES)


def entropy_bits(q):
    return -(q * np.log2(q)).sum(-1)


def kl_bits(p, r):
    """Row-wise directed KL in bits for strictly positive rows of equal shape."""
    return (p * (np.log2(p) - np.log2(r))).sum(-1)


def blahut(mass, cost, beta, tol=1e-12, max_iter=200000):
    """Blahut iteration that also returns the test channel (bits, slope beta)."""
    kernel = np.exp2(-beta * (cost - cost.min(1, keepdims=True)))
    prior = np.full(cost.shape[1], 1 / cost.shape[1])
    for _ in range(max_iter):
        z = kernel @ prior
        gain = (mass / z) @ kernel
        if np.log2(gain.max()) <= tol:
            break
        prior = prior * gain
        prior /= prior.sum()
    channel = kernel * prior / (kernel @ prior)[:, None]
    return channel


def information_bits(mass, channel):
    joint = mass[:, None] * channel
    out = joint.sum(0)
    ratio = np.divide(channel, out[None, :], out=np.ones_like(channel), where=(channel > 0) & (out[None, :] > 0))
    return float((joint * np.log2(ratio)).sum())


# --------------------------------------------------------------------------
# Block 1: the exact loss budget of a link
# --------------------------------------------------------------------------
def gray(index):
    return index ^ (index >> 1)


def index_channel(bits, error, labels):
    """Detection kernel of a fixed-width index sent over independent BSC(error) uses."""
    distance = np.array([[bin(int(a) ^ int(b)).count("1") for b in labels] for a in labels])
    return error**distance * (1 - error)**(bits - distance)


def conditional_mean(joint, q):
    """Posterior mean per output; an output of zero probability keeps a placeholder row."""
    mass = joint.sum(0)
    safe = np.where(mass > 0, mass, 1.0)
    mean = (joint / safe).T @ q
    mean[mass <= 0] = q[0]
    return mean


def loss_budget(beta=8.0, amplitude=0.2):
    n_x, n_u, bits = 64, 32, 5
    theta = 2 * np.pi * np.arange(n_x) / n_x
    q_x = circle(theta, amplitude)                       # P(Y | X)
    p_x = np.full(n_x, 1 / n_x)
    u_of_x = np.arange(n_x) // 2                         # frozen reference representation
    p_u = np.full(n_u, 1 / n_u)
    q_u = 0.5 * (q_x[0::2] + q_x[1::2])                  # Q = P(Y | U)
    cost = posterior_kl_bits(q_u, q_u)                   # object alphabet: the realizable posteriors
    w = blahut(p_u, cost, beta)                          # P(J | U)
    rows = []
    for name, snr_db, labels in [("clean", None, np.arange(n_u)),
                                 ("natural, 6 dB", 6.0, np.arange(n_u)),
                                 ("Gray, 6 dB", 6.0, gray(np.arange(n_u))),
                                 ("Gray, 2 dB", 2.0, gray(np.arange(n_u)))]:
        error = 0.0 if snr_db is None else 0.5 * erfc(np.sqrt(10 ** (snr_db / 10)))
        k = np.eye(n_u) if snr_db is None else index_channel(bits, error, labels)
        joint_uj = p_u[:, None] * w
        joint_uv = joint_uj @ k
        m_j = conditional_mean(joint_uj, q_u)            # E[Q | J]
        m_v = conditional_mean(joint_uv, q_u)            # E[Q | V]
        h_y_x = float(p_x @ entropy_bits(q_x))
        h_y_u = float(p_u @ entropy_bits(q_u))
        comp = float((joint_uj * (entropy_bits(m_j)[None, :])).sum() - h_y_u)   # I(Y;U|J)
        total_uv = float((joint_uv * (entropy_bits(m_v)[None, :])).sum() - h_y_u)  # I(Y;U|V)
        p_v = joint_uv.sum(0)
        decoders = {
            "lookup": q_u,                               # object stored under the received index
            "projection": q_u[np.argmin(posterior_kl_bits(m_v, q_u), axis=1)],
            "centroid": m_v,                             # unrestricted predictive receiver
        }
        for dec_name, table in decoders.items():
            dec = float(p_v @ kl_bits(m_v, table))
            direct = 0.0                                 # E[-log2 B_V(Y)] by full enumeration
            for x in range(n_x):
                pv_given_x = w[u_of_x[x]] @ k
                direct += p_x[x] * float(pv_given_x @ (-(q_x[x][None, :] * np.log2(table)).sum(1)))
            terms = dict(bayes=h_y_x, representation=h_y_u - h_y_x, compression=comp,
                         channel=total_uv - comp, decoder=dec)
            rows.append(dict(link=name, receiver=dec_name, bit_error=error,
                             information_bits=information_bits(p_u, w), **terms,
                             sum_of_terms=sum(terms.values()), direct_log_loss=direct,
                             identity_error=abs(sum(terms.values()) - direct)))
    return rows


# --------------------------------------------------------------------------
# Block 2 and 3: dense-alphabet receivers at equal information
# --------------------------------------------------------------------------
def fisher_speed_and_curvature(theta, amplitude=0.2):
    """Squared Fisher speed g and mixture curvature kappa of the posterior circle."""
    q = circle(theta, amplitude)
    d1 = -amplitude * np.sin(theta[:, None] - PHASES)
    d2 = -amplitude * np.cos(theta[:, None] - PHASES)
    g = (d1 * d1 / q).sum(1)
    along = (d2 * d1 / q).sum(1) / g
    normal = d2 - along[:, None] * d1
    kappa = np.sqrt((normal * normal / q).sum(1)) / g
    return g, kappa


def free_predictor(mass, q, beta, tol=1e-11, max_iter=20000):
    """Predictive receiver with unrestricted simplex outputs (free-centroid Blahut iteration).

    One movable reproduction per source state is initialised on the source itself.
    The iteration alternates the Gibbs rows, the output masses and the conditional
    means, so every limit point is a stationary point of I + beta * E KL(Q || b).
    It returns an achievable (information, distortion) pair in bits and the channel.
    """
    b = q.copy()
    prior = np.full(len(q), 1 / len(q))
    previous = np.inf
    for _ in range(max_iter):
        cost = posterior_kl_bits(q, np.maximum(b, 1e-300))
        kernel = np.exp2(-beta * (cost - cost.min(1, keepdims=True)))
        channel = kernel * prior
        channel /= channel.sum(1, keepdims=True)
        joint = mass[:, None] * channel
        prior = joint.sum(0)
        live = prior > 1e-300
        b[live] = (joint[:, live] / prior[live]).T @ q
        value = information_bits(mass, channel) + beta * float((joint * cost).sum())
        if previous - value < tol:
            break
        previous = value
    cost = posterior_kl_bits(q, np.maximum(b, 1e-300))
    joint = mass[:, None] * channel
    return information_bits(mass, channel), float((joint * cost).sum()), channel, b


def object_distortion_at_rate(mass, cost, rate, lo=1e-3, hi=1e5):
    """Distortion of the object receiver at a prescribed information, by bisection on the slope."""
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        point = finite_ba(mass, cost, mid, tol=1e-10, max_iter=20000)
        if point.rate < rate:
            lo = mid
        else:
            hi = mid
    a = finite_ba(mass, cost, lo, tol=1e-10, max_iter=20000)
    b = finite_ba(mass, cost, hi, tol=1e-10, max_iter=20000)
    if b.rate - a.rate < 1e-12:
        return 0.5 * (a.distortion + b.distortion)
    return a.distortion + (rate - a.rate) * (b.distortion - a.distortion) / (b.rate - a.rate)


def dense_curves(n=720, amplitude=0.2):
    theta = 2 * np.pi * np.arange(n) / n
    q = circle(theta, amplitude)
    p = np.full(n, 1 / n)
    eta, radius = 1 / 3 - amplitude, np.sqrt(1.5) * amplitude
    gamma_nat = (eta / (2 * radius)) ** 2                # curvature coefficient, nats
    cost_obj = posterior_kl_bits(q, q)
    g, kappa = fisher_speed_and_curvature(theta, amplitude)
    h_fisher_nat = float(np.log(2 * np.pi) + 0.5 * np.mean(np.log(g)))
    slopes = [12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048]
    obj, rows = [], []
    for beta in slopes:
        a = finite_ba(p, cost_obj, beta, tol=1e-10, max_iter=20000)
        obj.append((a.rate, a.distortion))
        rate, d_pred, w, b = free_predictor(p, q, beta)
        joint = p[:, None] * w
        mass = joint.sum(0)
        keep = mass > 1e-14
        m = b[keep]                                                    # centroids
        offset = m - CENTER
        norm = np.linalg.norm(offset, axis=1)
        projected = CENTER + radius * offset / norm[:, None]           # radial projection
        d_j = (joint[:, keep] * posterior_kl_bits(q, m)).sum(0) / mass[keep] * LN2   # nats
        proj_cost = kl_bits(m, projected) * LN2                                      # nats
        mean_d = float(mass[keep] @ d_j)
        d_obj = object_distortion_at_rate(p, cost_obj, rate)
        d_pred_nat = d_pred * LN2
        rows.append(dict(
            beta=beta, rate_bits=rate, D_predictive_bits=d_pred, D_object_bits=d_obj,
            gap_bits=d_obj - d_pred,
            gap_over_D2_nat=(d_obj - d_pred) * LN2 / d_pred_nat**2,
            curvature_lower_bits=gamma_nat * d_pred_nat**2 / LN2,
            local_law_bits=0.5 * float(np.mean(kappa**2)) * d_pred_nat**2 / LN2,
            projected_object_bits=(mean_d + float(mass[keep] @ proj_cost)) / LN2,
            projected_gap_over_D2_nat=float(mass[keep] @ proj_cost) / mean_d**2,
            reach_bound_bits=(mean_d + float(mass[keep] @ d_j**2) / (eta * radius**2)) / LN2,
            second_moment_ratio=float(mass[keep] @ d_j**2) / mean_d**2,
            asymptote_rate_bits=(h_fisher_nat - 0.5 * np.log(4 * np.pi * np.e * d_obj * LN2)) / LN2,
            asymptote_minus_rate_bits=(h_fisher_nat - 0.5 * np.log(4 * np.pi * np.e * d_obj * LN2)) / LN2 - rate))
    summary = dict(grid=n, amplitude=amplitude, eta=eta, radius=radius, gamma_nat=gamma_nat,
                   reach_coefficient_nat=1 / (eta * radius**2), h_fisher_nat=h_fisher_nat,
                   mean_kappa_sq=float(np.mean(kappa**2)), min_kappa_sq=float(kappa.min()**2),
                   max_kappa_sq=float(kappa.max()**2), half_mean_kappa_sq=0.5 * float(np.mean(kappa**2)))
    curves = dict(object=obj, predictive=[(r["rate_bits"], r["D_predictive_bits"]) for r in rows])
    return summary, rows, curves


def local_curvature_check(amplitude=0.2):
    """Projection cost over squared cell distortion for shrinking two-point cells."""
    rows = []
    for theta0 in (0.0, 0.4, np.pi / 3, 2.0):
        _, kappa = fisher_speed_and_curvature(np.array([theta0]), amplitude)
        for width in (0.4, 0.2, 0.1, 0.05):
            pts = circle(np.array([theta0 - width, theta0 + width]), amplitude)
            m = pts.mean(0, keepdims=True)
            d = float(np.mean(kl_bits(pts, np.repeat(m, 2, 0)))) * LN2
            fine = circle(theta0 + np.linspace(-0.2, 0.2, 400001), amplitude)
            cost = float(np.min(kl_bits(np.repeat(m, len(fine), 0), fine))) * LN2
            rows.append(dict(theta=float(theta0), half_width=width, ratio=cost / d**2,
                             half_kappa_sq=0.5 * float(kappa[0] ** 2)))
    return rows


# --------------------------------------------------------------------------
# Block 4: distortion floors at channel capacity, and the knowledge saving
# --------------------------------------------------------------------------
def opta(amplitude=0.2, n=24):
    theta = 2 * np.pi * np.arange(n) / n
    q = circle(theta, amplitude)
    p = np.full(n, 1 / n)
    cost = posterior_kl_bits(q, q)
    slopes = [0.25, 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    points = [finite_ba(p, cost, beta, tol=1e-9, max_iter=50000) for beta in slopes]
    rows = []
    for gain in (0.25, 0.5, 1.0, 2.0):
        capacity = float(np.log2(1 + gain**2 * 4.0 / 2))
        floor = max(0.0, max((a.dual - capacity) / a.beta for a in points))
        achievable = [a.distortion for a in points if a.rate <= capacity]
        rows.append(dict(gain=gain, capacity_bits=capacity, object_floor_bits=floor,
                         separation_achievable_bits=min(achievable) if achievable else float("nan")))
    tilted = [(1 + sign * 0.8 * np.cos(theta)) / n for sign in (-1, 1)]
    h_u = np.log2(n)
    h_u_k = float(np.mean([-(t * np.log2(t)).sum() for t in tilted]))
    saving = []
    for beta in (8, 16, 32, 64):
        plain = finite_ba(p, cost, beta, tol=1e-9, max_iter=50000)
        shared = conditional_ba([0.5, 0.5], tilted, [cost, cost], beta, tol=1e-9, max_iter=50000)
        saving.append(dict(beta=beta, D_plain=plain.distortion, D_shared=shared.distortion,
                           rate_plain=plain.rate, rate_shared=shared.rate,
                           lagrangian_saving_bits=plain.rate + beta * plain.distortion
                           - shared.rate - beta * shared.distortion))
    return rows, dict(mutual_information_U_K_bits=h_u - h_u_k, slopes=saving)


def outage_example(rate_bits, uses):
    """Rayleigh outage probability of a quasi-static complex link for a required source rate."""
    rows = []
    for snr_db in (0, 5, 10, 15, 20, 25, 30):
        snr = 10 ** (snr_db / 10)
        rows.append(dict(snr_db=snr_db, outage=float(1 - np.exp(-(2 ** (rate_bits / uses) - 1) / snr))))
    return rows


def plot(summary, rows, curves, path):
    import matplotlib.pyplot as plt

    obj, pred = np.array(curves["object"]), np.array(curves["predictive"])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.7), layout="constrained")
    d = np.geomspace(2e-4, 0.2, 200)
    asym = (summary["h_fisher_nat"] - 0.5 * np.log(4 * np.pi * np.e * d * LN2)) / LN2
    ax.semilogx(obj[:, 1], obj[:, 0], color="#124971", label="Object receiver")
    ax.semilogx(pred[:, 1], pred[:, 0], "--", color="#a32437", label="Predictive receiver")
    ax.semilogx(d, asym, ":", color=".3", label=r"$h_F-\frac{1}{2}\ln(4\pi eD)$")
    ax.set(xlabel="KL distortion $D$ (bits)", ylabel="Information (bits/state)", xlim=(2e-4, 0.2),
           ylim=(0, 6.5), title="(a) Task-rank law, $r=1$")
    ax.legend(frameon=False, fontsize=7)
    dp = np.array([r["D_predictive_bits"] for r in rows])
    gap = np.array([r["gap_bits"] for r in rows])
    bx.loglog(dp, gap, "o-", color="#124971", ms=3, label="Computed gap")
    bx.loglog(dp, [r["local_law_bits"] for r in rows], ":", color=".3",
              label=r"$\frac{1}{2}\overline{\kappa_m^2}\,D^2$")
    bx.loglog(dp, [r["curvature_lower_bits"] for r in rows], "--", color="#a32437", label="Lower bound, Theorem 2")
    bx.loglog(dp, [r["projected_object_bits"] - r["D_predictive_bits"] for r in rows], "-",
              color="#6a8f3c", lw=1, label="Projected centroids")
    bx.set(xlabel=r"$D_{\rm pred}(R)$ (bits)", ylabel=r"$D_{\rm obj}(R)-D_{\rm pred}(R)$ (bits)",
           title="(b) Price of an object")
    bx.legend(frameon=False, fontsize=7)
    for axis in (ax, bx):
        axis.grid(alpha=.2, which="both")
        axis.title.set_fontsize(9)
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
    parser.add_argument("--grid", type=int, default=720)
    args = parser.parse_args()
    budget = loss_budget()
    worst = max(r["identity_error"] for r in budget)
    if worst > 1e-10:
        raise ArithmeticError(f"loss budget identity violated by {worst}")
    summary, rows, curves = dense_curves(args.grid)
    floors, knowledge = opta()
    record = dict(loss_budget=budget, dense_summary=summary, equal_information=rows,
                  dense_curves=curves, local_curvature=local_curvature_check(),
                  capacity_floors=floors, knowledge=knowledge,
                  outage_example=outage_example(rate_bits=1.2601, uses=1))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "theory_checks.json").write_text(json.dumps(record, indent=1))
    if args.figure:
        plot(summary, rows, curves, args.figure)
    print(f"loss-budget identity error {worst:.2e}; wrote {args.output / 'theory_checks.json'}")


if __name__ == "__main__":
    main()
