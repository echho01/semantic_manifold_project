"""Self-test of tda.py against persim / ripser and known closed forms.

Run:  .venv/bin/python code/topology/selftest_tda.py
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np

import tda



def run_checks():
    rng = np.random.default_rng(0)
    ok = True
    metrics = {}

    # 1. bottleneck vs persim on random finite diagrams
    import persim
    err = 0.0
    for _ in range(300):
        dg = []
        for _ in range(2):
            m = int(rng.integers(1, 15))
            a = rng.uniform(0, 1, (m, 2))
            a[:, 1] = a[:, 0] + rng.exponential(0.3, m)
            dg.append(a)
        err = max(err, abs(tda.bottleneck(*dg) - persim.bottleneck(*dg)))
    print(f"bottleneck  max|mine - persim|            = {err:.2e}")
    ok &= err < 1e-8
    metrics["bottleneck_vs_persim_max_abs_error"] = float(err)
    metrics["bottleneck_random_pairs"] = 300

    # 2. H0 deaths from the MST vs ripser (radius convention)
    X = rng.normal(size=(300, 5))
    d = tda.rips(X, maxdim=0)[0]
    d = np.sort(d[np.isfinite(d[:, 1]), 1])[::-1]
    e = np.abs(d - tda.h0_deaths(X)).max()
    print(f"H0 deaths   max|MST - ripser|             = {e:.2e}")
    ok &= e < 1e-6
    metrics["h0_mst_vs_ripser_max_abs_error"] = float(e)

    # 3. Cech H1 of an ellipse (a, b) = (1, 0.5): bar [~0, b), reach b^2/a
    th = np.linspace(0, 2 * np.pi, 1500, endpoint=False)
    A = tda.alpha2d(np.c_[np.cos(th), 0.5 * np.sin(th)])
    h1 = A[1][np.argmax(tda.persistence(A[1]))]
    print(f"ellipse     Cech H1 bar = [{h1[0]:.4f}, {h1[1]:.4f})   "
          f"predicted [0, 0.5); reach 0.25")
    ok &= abs(h1[1] - 0.5) < 1e-3 and h1[0] < 5e-3
    metrics["ellipse_alpha_h1_interval"] = h1.tolist()

    # 4. alpha (Cech) H0 equals Rips H0
    P = rng.uniform(size=(400, 2))
    A = tda.alpha2d(P)
    h0 = np.sort(A[0][np.isfinite(A[0][:, 1]), 1])[::-1]
    e = np.abs(h0 - tda.h0_deaths(P)).max()
    print(f"alpha H0    max|alpha - MST|              = {e:.2e}")
    ok &= e < 1e-9
    metrics["alpha_h0_vs_mst_max_abs_error"] = float(e)

    # 5. alpha H1 vs ripser H1 on a planar circle: Cech death = rho,
    #    Rips death (radius) = rho*sqrt(3)/2 (equilateral triangle)
    pts = np.c_[np.cos(th[::5]), np.sin(th[::5])]
    c = tda.alpha2d(pts)[1]
    r = tda.rips(pts, maxdim=1)[1]
    cd = c[np.argmax(tda.persistence(c)), 1]
    rd = r[np.argmax(tda.persistence(r)), 1]
    print(f"circle      Cech death {cd:.4f} (1.0),  Rips death {rd:.4f} "
          f"({np.sqrt(3) / 2:.4f})")
    ok &= abs(cd - 1) < 1e-3 and abs(rd - np.sqrt(3) / 2) < 2e-2
    metrics["circle_alpha_death"] = float(cd)
    metrics["circle_rips_death"] = float(rd)
    native_error = 0.0
    for n in [3, 8, 25, 60]:
        for _ in range(5):
            X = rng.normal(size=(n, 2))
            native = tda.alpha2d(X, backend="native")
            independent = tda.alpha2d(X, backend="gudhi")
            native_error = max(native_error, max(tda.bottleneck_by_degree(native, independent)))
    ok &= native_error < 1e-10
    metrics["native_alpha_vs_gudhi_max_bottleneck_error"] = float(native_error)
    metrics["native_alpha_clouds"] = 20
    print(f"native alpha vs GUDHI max dB             = {native_error:.2e}")

    print("ALL PASS" if ok else "FAILURES ABOVE")
    return {"passed": bool(ok), "metrics": metrics,
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), Path(tda.__file__)]},
        "dependencies": {name: importlib.metadata.version(name) for name in ["numpy", "scipy", "ripser", "persim", "gudhi"]}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    report = run_checks()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
