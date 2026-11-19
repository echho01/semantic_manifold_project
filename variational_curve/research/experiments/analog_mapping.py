"""Analog mapping of a nonuniform circular source: area, reach and the companding law.

A source angle theta with density p on the circle, semantic metric g = dtheta^2, is sent as
F(theta) = rho (cos s(theta), sin s(theta)) over two real Gaussian uses with noise variance
sigma^2. The receiver returns the nearest point of the surface. Every such surface is the same
circle of radius rho, so power (rho^2), reach (rho) and length (2 pi rho) are identical and only
the speed s'(theta), the area density of the surface, differs:

  isometric   s' constant;
  companded   s' proportional to p^(1/3), the equality case of the area bound for m = 1.

The semantic distortion E (theta - theta_hat)^2 is computed exactly by quadrature of the angular
error density of a two-dimensional Gaussian. It is compared with the weak-noise formula
sigma^2 E[1/|F'|^2] and with the bound sigma^2 (int p^(1/3))^3 / L^2 of the theorem.

    python experiments/analog_mapping.py --output results/analog_mapping.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import erfc


def angular_error_density(phi, snr):
    """Density of the phase of (sqrt(snr), 0) + N(0, I/2 * 2) scaled so that rho/sigma = sqrt(snr)."""
    a = np.sqrt(snr)                                   # rho / sigma
    c, s2 = np.cos(phi), np.sin(phi) ** 2
    return (np.exp(-a * a / 2) / (2 * np.pi)
            + a * c * np.exp(-a * a * s2 / 2) * (1 - 0.5 * erfc(a * c / np.sqrt(2))) / np.sqrt(2 * np.pi))


def surface(p, exponent, n):
    """Reparametrization s(theta) with s' proportional to p^exponent and total turn 2 pi."""
    theta = (np.arange(n) + 0.5) * 2 * np.pi / n
    speed = p(theta) ** exponent
    speed *= n / speed.sum()                           # mean angular speed one
    edges = np.concatenate(([0.0], np.cumsum(speed))) * 2 * np.pi / n
    return theta, speed, edges


def exact_distortion(p, exponent, snr, n=4000, n_phi=4001):
    theta, speed, edges = surface(p, exponent, n)
    s_mid = 0.5 * (edges[1:] + edges[:-1])
    grid_theta = np.arange(n + 1) * 2 * np.pi / n      # theta at the cell edges
    phi = np.linspace(-np.pi, np.pi, n_phi)
    f = angular_error_density(phi, snr)
    w = np.full(n_phi, phi[1] - phi[0])
    w[[0, -1]] /= 2
    f = f / np.sum(f * w)
    weights = p(theta) / p(theta).sum()
    total = 0.0
    for i in range(0, n, 8):                           # Subsample source points; the integrand is smooth.
        received = np.mod(s_mid[i] + phi, 2 * np.pi)
        decoded = np.interp(received, edges, grid_theta)
        err = np.abs(decoded - theta[i])
        err = np.minimum(err, 2 * np.pi - err)
        total += weights[i] * 8 * np.sum(f * w * err**2)
    return float(total)


def run(tilt=0.85):
    p = lambda t: (1 + tilt * np.cos(t)) / (2 * np.pi)
    theta = (np.arange(200000) + 0.5) * 2 * np.pi / 200000
    d = theta[1] - theta[0] + 0 * theta
    norm = float(np.sum(p(theta) ** (1 / 3) * d) ** 3)           # ||p||_{1/3}
    rows = []
    for snr_db in (10, 15, 20, 25, 30, 40):
        snr = 10 ** (snr_db / 10)                                # rho^2 / sigma^2
        row = dict(snr_db=snr_db)
        for name, exponent in (("isometric", 0.0), ("companded", 1 / 3)):
            _, speed, _ = surface(p, exponent, 200000)
            # Weak-noise formula sigma^2 E[1/|F'|^2] with |F'| = rho s', in units of radian^2.
            weak = float(np.sum(p(theta) / speed**2 * d) / snr)
            row[name] = dict(exact=exact_distortion(p, exponent, snr), weak_noise=weak)
        row["area_bound"] = norm / (2 * np.pi) ** 2 / snr        # sigma^2 ||p||_{1/3} / L^2
        row["exact_ratio"] = row["companded"]["exact"] / row["isometric"]["exact"]
        rows.append(row)
    return dict(tilt=tilt, norm_one_third=norm, predicted_ratio=norm / (2 * np.pi) ** 2,
                tube_bound_over_length=2.25, rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "results" / "analog_mapping.json")
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=1) + "\n")
    print("predicted ratio", round(result["predicted_ratio"], 4))
    for r in result["rows"]:
        print(r["snr_db"], "iso %.3e (weak %.3e)" % (r["isometric"]["exact"], r["isometric"]["weak_noise"]),
              "comp %.3e (weak %.3e)" % (r["companded"]["exact"], r["companded"]["weak_noise"]),
              "bound %.3e" % r["area_bound"], "ratio %.4f" % r["exact_ratio"])
