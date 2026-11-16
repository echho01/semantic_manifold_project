"""Intrinsic-dimension and local Monge-chart geometry diagnostics.

A local curvature radius is not a lower bound on global reach: a distant
sheet may be arbitrarily close. Full distance matrices allocate n*m output.
"""
from __future__ import annotations
import numpy as np
from sklearn.neighbors import NearestNeighbors


def _cloud(X):
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2 or len(X) < 3 or not np.isfinite(X).all():
        raise ValueError("expected a finite [n>=3, ambient_dimension] cloud")
    return X


def pair_distances(X, Y=None, block=1024):
    """Squared Euclidean matrix; temporaries blocked, output remains n*m."""
    X = _cloud(X)
    Y = X if Y is None else np.asarray(Y, dtype=np.float64)
    out = np.empty((len(X), len(Y)))
    yn = (Y * Y).sum(1)
    for i in range(0, len(X), block):
        a = X[i:i + block]
        out[i:i + block] = np.maximum((a * a).sum(1)[:, None] + yn - 2 * a @ Y.T, 0)
    return out


def knn(X, k, block=4096):
    """Exact k nearest neighbours, excluding the query's own row index."""
    X = _cloud(X)
    if not 1 <= k < len(X):
        raise ValueError("require 1 <= k < number of samples")
    distances, indices = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    out_i = np.empty((len(X), k), dtype=int)
    out_d = np.empty((len(X), k))
    for i, (ii, dd) in enumerate(zip(indices, distances)):
        mask = ii != i
        out_i[i] = ii[mask][:k]
        out_d[i] = dd[mask][:k] ** 2
    return out_i, out_d


def twonn_dimension(X, frac=.9, max_points=20000, seed=0):
    """TwoNN through-origin fit using the full-sample empirical CDF."""
    X = _cloud(X)
    if not 0 < frac <= 1 or max_points < 20:
        raise ValueError("invalid TwoNN settings")
    original = len(X)
    X = np.unique(X, axis=0)
    unique = len(X)
    if unique < 20:
        raise ValueError("TwoNN requires >=20 distinct embeddings; class prompts are discrete prototypes")
    rng = np.random.default_rng(seed)
    if len(X) > max_points:
        X = X[rng.choice(len(X), max_points, replace=False)]
    _, d2 = knn(X, 2)
    valid = (d2[:, 0] > 0) & np.isfinite(d2).all(1)
    mu = np.sort(np.sqrt(d2[valid, 1] / d2[valid, 0]))
    n = len(mu)
    m = max(2, min(n, int(frac * n)))
    x = np.log(mu[:m])
    y = -np.log1p(-(np.arange(1, m + 1) - .5) / n)
    if x @ x <= 0:
        raise ValueError("nearest-neighbour ratios do not identify a dimension")
    d = float((x @ y) / (x @ x))
    return {"dim": d, "n": m, "n_total": n, "n_original": original,
            "n_unique": unique, "duplicates_removed": original - unique,
            "residual": float(np.linalg.norm(y - d * x) / np.linalg.norm(y))}


def local_curvature(X, pca_dim, k_nn=120, n_anchor=4000, seed=0):
    """Fit II with a symmetric quadratic basis and the induced chart metric.

    Normal displacements are fitted in ambient coordinates, avoiding a
    missing normal basis when ambient dimension exceeds neighbour count.
    The graph Hessian is projected to the fitted graph's normal space,
    then orthonormalized using its Jacobian metric.
    """
    X = _cloud(X)
    n, ambient = X.shape
    d = int(pca_dim)
    params = 1 + d + d * (d + 1) // 2
    if not 1 <= d < ambient:
        raise ValueError("pca_dim must be in [1, ambient_dimension-1]")
    if not params < k_nn < n:
        raise ValueError(f"k_nn must exceed {params} quadratic parameters and be smaller than n={n}")
    anchors = np.random.default_rng(seed).choice(n, min(n_anchor, n), replace=False)
    _, near = NearestNeighbors(n_neighbors=k_nn + 1).fit(X).kneighbors(X[anchors])
    scalar, ii2, h2, radii, fit_errors = [], [], [], [], []
    pairs = [(i, j) for i in range(d) for j in range(i, d)]
    for anchor, ids in zip(anchors, near):
        ids = ids[ids != anchor][:k_nn]
        v = X[ids] - X[anchor]
        _, _, vt = np.linalg.svd(v - v.mean(0), full_matrices=False)
        e = vt[:d].T
        xi = v @ e
        rho = v - xi @ e.T
        scale = np.sqrt(np.mean(np.sum(xi * xi, axis=1)))
        if scale <= 0:
            raise ValueError("degenerate local neighbourhood")
        z = xi / scale
        quadratic = np.column_stack([.5 * z[:, i] ** 2 if i == j else z[:, i] * z[:, j] for i, j in pairs])
        design = np.column_stack([np.ones(len(z)), z, quadratic])
        coef, _, rank, _ = np.linalg.lstsq(design, rho, rcond=1e-10)
        if rank < params:
            raise ValueError("rank-deficient quadratic fit; lower pca_dim or increase neighbours")
        jac = e + coef[1:1 + d].T / scale
        val, vec = np.linalg.eigh(jac.T @ jac)
        inv = (vec / val) @ vec.T
        inverse_root = (vec / np.sqrt(val)) @ vec.T
        hessian = np.zeros((d, d, ambient))
        for coefficient, (i, j) in zip(coef[1 + d:], pairs):
            hessian[i, j] = hessian[j, i] = coefficient / (scale * scale)
        flat = hessian.reshape(-1, ambient)
        normal = flat - (flat @ jac) @ inv @ jac.T
        h = np.einsum("ia,jb,ijc->abc", inverse_root, inverse_root, normal.reshape(d, d, ambient))
        ii = float(np.sum(h * h))
        trace = np.einsum("iic->c", h)
        hs = float(trace @ trace)
        scalar.append(hs - ii)
        ii2.append(ii)
        h2.append(hs)
        radii.append(np.sqrt(d / ii) if ii > 1e-24 else np.inf)
        fit_errors.append(float(np.linalg.norm(rho - design @ coef) / max(np.linalg.norm(rho), 1e-15)))
    scalar, ii2, h2, radii = map(np.asarray, (scalar, ii2, h2, radii))
    finite = radii[np.isfinite(radii)]
    return {"scal_mean": float(scalar.mean()), "scal_quantiles": np.quantile(scalar, [.1, .5, .9]).tolist(),
            "ii2_mean": float(ii2.mean()), "h2_mean": float(h2.mean()),
            "ii2_quantiles": np.quantile(ii2, [.1, .5, .9]).tolist(),
            "curvature_radius_median": float(np.median(finite)) if len(finite) else None,
            "curvature_radius_interpretation": "sqrt(d/||II||^2), local RMS radius; not global reach",
            "fit_relative_residual_median": float(np.median(fit_errors)),
            "d": d, "k_nn": k_nn, "n_anchor": len(anchors),
            "anchor_indices": anchors.tolist(), "scalar_local": scalar.tolist(),
            "ii2_local": ii2.tolist(), "h2_local": h2.tolist()}


def local_scalar_profile(X, pca_dim, k_nn=64, n_anchor=2000, seed=0):
    return local_curvature(X, pca_dim, k_nn, n_anchor, seed)
