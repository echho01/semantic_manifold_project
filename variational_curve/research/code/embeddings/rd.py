"""Held-out fixed-rate vector quantization of embedding clouds.

The rate is log2(codebook size). Its second-order residual is an empirical
diagnostic, not a Shannon RD coefficient. Continuous RD coefficients are
validated by code/synthetic/run_paper1.py.
"""
from __future__ import annotations
import numpy as np
from sklearn.cluster import MiniBatchKMeans


def split_indices(n, seed=0, test_fraction=.2):
    if n < 10 or not 0 < test_fraction < 1:
        raise ValueError("need n>=10 and test_fraction in (0,1)")
    order = np.random.default_rng(seed).permutation(n)
    nt = max(2, int(round(n * test_fraction)))
    return order[nt:], order[:nt]


def quantize_curve(X, sizes, n_fit=200000, batch=8192, iters=50, seed=0,
                   geodesic=False, tau=None, test_fraction=.2):
    """Train on one split and measure ambient MSE on another.

    A scalar curvature/reach proxy cannot supply geodesic ground truth for
    off-manifold k-means centroids. Actual geodesic evaluation requires a
    specified manifold metric and manifold-valued reproductions.
    """
    if geodesic:
        raise ValueError("geodesic evaluation requires a specified manifold metric; no reach-based surrogate is used")
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2 or not np.isfinite(X).all():
        raise ValueError("invalid embedding cloud")
    train_idx, test_idx = split_indices(len(X), seed, test_fraction)
    test_idx = test_idx[:n_fit]
    train, test = X[train_idx], X[test_idx]
    out = {"sizes": [], "D": [], "median_D": [], "D_standard_error": [], "skipped_sizes": [],
           "train_count": len(train), "test_count": len(test), "test_indices": test_idx.tolist(),
           "distortion": "ambient squared Euclidean", "rate": "fixed-rate log2(N) bits",
           "reproduction_alphabet": "unconstrained ambient centroids"}
    for requested in sizes:
        n = int(requested)
        if n < 1 or n > len(train):
            out["skipped_sizes"].append({"size": n, "reason": "outside available training-sample range"})
            continue
        km = MiniBatchKMeans(n_clusters=n, batch_size=max(batch, n), n_init=3,
                             max_iter=iters, random_state=seed).fit(train)
        labels = km.predict(test)
        d2 = np.sum((test - km.cluster_centers_[labels]) ** 2, axis=1)
        out["sizes"].append(n)
        out["D"].append(float(d2.mean()))
        out["median_D"].append(float(np.median(d2)))
        out["D_standard_error"].append(float(d2.std(ddof=1) / np.sqrt(len(d2))))
    return out


def fit_second_order(curve, d):
    n = np.asarray(curve["sizes"], float)
    dist = np.asarray(curve["D"], float)
    mask = (dist > 0) & np.isfinite(dist) & (n > 0)
    n, dist = n[mask], dist[mask]
    if len(dist) < 3:
        return {"status": "insufficient positive-distortion points"}
    y = np.log2(n) + (d / 2) * np.log2(dist)
    scale = dist.max()
    a = np.column_stack([np.ones(len(dist)), dist / scale, (dist / scale) ** 2])
    coef = np.linalg.lstsq(a, y, rcond=None)[0]
    res = y - a @ coef
    return {"intercept_bits": float(coef[0]), "slope_bits_per_D": float(coef[1] / scale),
            "quadratic_bits_per_D2": float(coef[2] / scale ** 2),
            "rms_bits": float(np.sqrt(np.mean(res * res))), "D": dist.tolist(), "residual_bits": y.tolist(),
            "interpretation": "empirical fixed-rate residual; does not estimate the Shannon RD curvature coefficient"}


def predict_coefficients(scal_mean, d, ii2_mean=None, h2_mean=None):
    """Continuous manifold-constrained Shannon RD coefficients in nats."""
    if d <= 0:
        raise ValueError("dimension must be positive")
    result = {"geodesic_nats_per_D": scal_mean / (6 * d)}
    if ii2_mean is not None:
        h2 = ii2_mean + scal_mean if h2_mean is None else h2_mean
        result["chordal_nats_per_D"] = (h2 - 2 * ii2_mean) / (8 * d)
    return result


def local_profile(X, labels, centers, scal_local, rng_seed=0):
    """Within-cell centered distortion/curvature correlation diagnostic."""
    from scipy.stats import spearmanr
    labels = np.asarray(labels)
    d2 = np.sum((X - centers[labels]) ** 2, axis=1)
    scal = np.asarray(scal_local, dtype=float)
    if len(scal) != len(d2):
        raise ValueError("curvature values must align with points")
    residual_d, residual_s = d2.copy(), scal.copy()
    for label in np.unique(labels):
        mask = labels == label
        residual_d[mask] -= d2[mask].mean()
        residual_s[mask] -= scal[mask].mean()
    rho, p = spearmanr(residual_d, residual_s)
    return {"spearman_rho": float(rho), "p_value": float(p), "interpretation": "within-cell descriptive correlation"}
