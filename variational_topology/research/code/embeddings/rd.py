"""Rate-distortion experiments on real embeddings.

Theory under test (Paper 1, Theorem 3.1 and Proposition 3.4):

  geodesic criterion:   R(D) = h - (d/2) log2(2 pi e D / d)
                               + Sbar / (6 d ln 2) D + o(D)
  ambient criterion:    R(D) = h - (d/2) log2(2 pi e D / d)
                               - (1/16) E[(k1-k2)^2] / ln 2 D + o(D)

The rate here is log2 N bits per codebook; the entropy h is a constant
that does not affect the slope, so the test is on the coefficient of D:
write D(N) = c N^{-2/d} (1 + b N^{-2/d}), so

    rel(N) := log2(D(N) N^{2/d} / c) = b N^{-2/d} + O(N^{-4/d}),

and b is the curvature coefficient in nats:  b = +Scal/(6d) for the
geodesic criterion, b = -(1/16) E[(k1-k2)^2] for the ambient one.

Two estimators are produced:
  (i)  aggregate: fit rel(N) = a + b/n over the codebook ladder, where
       n = N^{2/d}.  `a` absorbs the lattice-class constant, `b` is the
       curvature term.  This is the same regression used on synthetic
       manifolds, so the numbers are directly comparable.
  (ii) local:    at fixed N, correlate local distortion (excess over the
       cell mean) with the local scalar curvature estimate.  The
       prediction is a signed linear relation.
"""
from __future__ import annotations

import numpy as np
from sklearn.cluster import MiniBatchKMeans

from geometry import pair_distances


def quantize_curve(X: np.ndarray, sizes: list[int], n_fit: int = 200000,
                   batch: int = 8192, iters: int = 50, seed: int = 0,
                   geodesic: bool = True, tau: float | None = None) -> dict:
    """Train k-means codebooks and measure mean squared distortion.

    For `geodesic`, the chordal squared distance is converted with the
    chord-arc expansion of Lemma 2.2,
        d_g^2 = d_E^2 (1 + d_E^2 / (3 tau^2)),
    which is the regime the theory describes (D << tau^2).
    """
    rng = np.random.default_rng(seed)
    n = len(X)
    fit_idx = rng.choice(n, min(n_fit, n), replace=False)
    Xf = X[fit_idx]
    out = {"sizes": [], "D": [], "within": [], "labels": []}
    for N in sizes:
        km = MiniBatchKMeans(n_clusters=N, batch_size=batch, n_init=3,
                             max_iter=iters, random_state=seed).fit(X)
        lab = km.predict(Xf)
        d2 = ((Xf - km.cluster_centers_[lab]) ** 2).sum(1)
        if geodesic and tau is not None:
            d2 = d2 * (1.0 + d2 / (3.0 * tau ** 2))
        out["sizes"].append(int(N))
        out["D"].append(float(d2.mean()))
        out["within"].append(float(np.median(d2)))
        out["labels"].append(lab)
    out["fit_idx"] = fit_idx
    return out


def fit_second_order(curve: dict, d: float) -> dict:
    """Fit the residual of the rate--distortion curve.

    For the entropy-coded optimum, R(D) + (d/2) log2 D is constant plus a
    term linear in D.  For a fixed-rate (k-means) curve, the leading
    constant differs from log2(2 pi e / d) by the lattice-class constant,
    which is unknown for a learned embedding.  The intercept of the fit
    absorbs it; the slope is the curvature coefficient.

    Model:  y(N) := log2 N + (d/2) log2 D(N) = a + b D(N) + c D(N)^2.
    """
    N = np.asarray(curve["sizes"], float)
    D = np.asarray(curve["D"], float)
    y = np.log2(N) + (d / 2.0) * np.log2(D)
    A = np.vstack([np.ones_like(D), D, D ** 2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b, c = (float(v) for v in coef)
    resid = y - A @ coef
    return {"a": a, "b": b, "c": c,
            "rms": float(np.sqrt((resid ** 2).mean())),
            "n": N.tolist(), "D": D.tolist(), "y": y.tolist()}


def predict_coefficients(scal_mean: float, d: float,
                         anisotropy: float | None = None) -> dict:
    """Closed-form predictions in nats, matching Theorem 3.1 / Prop. 3.4."""
    out = {"geodesic": scal_mean / (6.0 * d)}
    if anisotropy is not None:
        out["chordal"] = -0.25 * anisotropy
    return out


def local_profile(X: np.ndarray, labels: np.ndarray, centers: np.ndarray,
                  scal_local: np.ndarray, rng_seed: int = 0) -> dict:
    """Correlate local excess distortion with local scalar curvature.

    `scal_local` is the per-point curvature estimate (one value per row of
    X, e.g. broadcast from a local fit or an interpolated field).  The
    prediction is a positive correlation under the geodesic criterion when
    the excess is measured against the cell mean at the same cell scale.
    """
    from scipy.stats import spearmanr
    d2 = ((X - centers[labels]) ** 2).sum(1)
    rho, p = spearmanr(d2, scal_local)
    return {"spearman_rho": float(rho), "p_value": float(p)}
