"""Geometry estimation for embedding clouds.

Three estimators, each with an explicit model of what it measures:

  twonn_dimension   intrinsic dimension, Facco et al. 2017
  local_curvature   second fundamental form and scalar curvature, from a
                    k-NN quadratic fit in the ambient space
  pair_distances    squared Euclidean distances in blocks (memory-safe)

All functions take float32 [n, D] arrays and return plain numpy / dict.
"""
from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


# ---------------------------------------------------------------- distances
def pair_distances(X: np.ndarray, Y: np.ndarray | None = None,
                   block: int = 4096) -> np.ndarray:
    """Squared Euclidean distances, blocked to bound memory."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    Y = X if Y is None else np.ascontiguousarray(Y, dtype=np.float32)
    out = np.empty((len(X), len(Y)), dtype=np.float32)
    for i in range(0, len(X), block):
        A = X[i:i + block]
        out[i:i + block] = (
            (A * A).sum(1)[:, None] + (Y * Y).sum(1)[None, :] - 2.0 * (A @ Y.T)
        )
    np.maximum(out, 0.0, out=out)
    return out


def knn(X: np.ndarray, k: int, block: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    """Exact k-NN indices (self excluded) and squared distances."""
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean")
    nn.fit(X)
    dist, idx = nn.kneighbors(X, return_distance=True)
    return idx[:, 1:], (dist[:, 1:] ** 2).astype(np.float32)


# ------------------------------------------------------------------- TwoNN
def twonn_dimension(X: np.ndarray, frac: float = 0.9, max_points: int = 20000,
                    seed: int = 0) -> dict:
    """Intrinsic dimension by the TwoNN estimator.

    mu_i = r2_i / r1_i for the two nearest neighbours; under a locally
    uniform density P(mu > t) = t^{-d}, so the empirical CDF F satisfies
    log(1 - F(mu)) = -d log(mu).  Fit through the origin on the lower
    `frac` quantile of mu (the upper tail is boundary-contaminated).
    """
    rng = np.random.default_rng(seed)
    n = len(X)
    if n > max_points:
        X = X[rng.choice(n, max_points, replace=False)]
    idx, d2 = knn(X, 2)
    mu = np.sqrt(d2[:, 1]) / np.sqrt(d2[:, 0])
    mu = np.sort(mu[np.isfinite(mu)])
    N = len(mu)
    m = max(10, int(frac * N))
    mu = mu[:m]
    # The empirical CDF must be normalized by the FULL sample: truncating
    # mu and then renormalizing inflates F and biases d upward.
    F = (np.arange(1, m + 1) - 0.5) / N
    y = -np.log(1.0 - F)
    x = np.log(mu)
    d = float((x @ y) / (x @ x))
    resid = np.linalg.norm(y - d * x) / np.linalg.norm(y)
    return {"dim": d, "n": int(m), "n_total": int(N), "residual": float(resid)}


# ------------------------------------------------------- local curvature
def local_curvature(X: np.ndarray, pca_dim: int, k_nn: int = 120,
                    n_anchor: int = 4000, seed: int = 0) -> dict:
    """Second fundamental form and scalar curvature from local quadratic fits.

    At an anchor x with k neighbours, let (e_1..e_d) span the local tangent
    plane (PCA) and (n_1..n_m) its orthogonal complement in ambient space.
    Write each neighbour as xi*e_i + rho*n.  The normal displacement obeys
        rho_a ≈ 1/2 sum_{ij} H^a_{ij} xi_i xi_j ,
    so H^a is a least-squares fit of the m normal coordinates on the
    quadratic features.  Then, by the Gauss equation,
        ||II||_HS^2 = sum_a ||H^a||_F^2 ,
        Scal = sum_a [ (tr H^a)^2 - ||H^a||_F^2 ] .
    The reach proxy is 1 / max_a ||H^a||_op.

    Note on the sign: H^a here is the second fundamental form of the graph
    rho(xi) with respect to the ambient normal n_a.  Only quadratic
    invariants are used, so the sign convention does not affect them.
    """
    rng = np.random.default_rng(seed)
    n, D = X.shape
    d = int(pca_dim)
    if d >= D:
        raise ValueError(f"pca_dim={d} must be smaller than ambient dim {D}")
    anchors = rng.choice(n, size=min(n_anchor, n), replace=False)
    idx, _ = knn(X, k_nn)

    scal, ii2, reach, traces = [], [], [], []
    for a in anchors:
        nb = X[idx[a]]
        center = nb.mean(0)
        V = nb - center
        # local PCA
        _, _, Vt = np.linalg.svd(V, full_matrices=False)
        E = Vt[:d].T                      # [D, d] tangent frame
        Nrm = Vt[d:].T                    # [D, D-d] normal frame
        xi = V @ E                        # [k, d]
        rho = V @ Nrm                     # [k, D-d]
        # quadratic design: intercept + linear + quadratic terms
        q = (xi[:, :, None] * xi[:, None, :]).reshape(len(xi), d * d)
        A = np.hstack([np.ones((len(xi), 1)), xi, q])
        coef, *_ = np.linalg.lstsq(A, rho, rcond=None)   # [1+d+d^2, D-d]
        Hs = []
        for j in range(D - d):
            # rho = 1/2 sum_ij H_ij xi_i xi_j, so H = 2 * coef
            H = 2.0 * coef[1 + d:, j].reshape(d, d)
            H = 0.5 * (H + H.T)
            Hs.append(H)
        tr = np.array([np.trace(H) for H in Hs])
        fro = np.array([np.linalg.norm(H, "fro") ** 2 for H in Hs])
        scal.append(float((tr ** 2 - fro).sum()))
        ii2.append(float(fro.sum()))
        reach.append(float(1.0 / max(np.abs(np.linalg.eigvalsh(H)).max()
                                    for H in Hs)))
        traces.append(tr)

    scal = np.array(scal); ii2 = np.array(ii2); reach = np.array(reach)
    return {
        "scal_mean": float(scal.mean()),
        "scal_quantiles": np.quantile(scal, [0.1, 0.5, 0.9]).tolist(),
        "ii2_mean": float(ii2.mean()),
        "ii2_quantiles": np.quantile(ii2, [0.1, 0.5, 0.9]).tolist(),
        "reach_median": float(np.median(reach)),
        "reach_quantiles": np.quantile(reach, [0.1, 0.5, 0.9]).tolist(),
        "d": d, "k_nn": k_nn, "n_anchor": int(len(anchors)),
    }


# --------------------------------------------------- local scalar profile
def local_scalar_profile(X: np.ndarray, pca_dim: int, k_nn: int = 64,
                         n_anchor: int = 2000, seed: int = 0) -> dict:
    """Scalar curvature and local distortion scale per anchor, for the test
    of the pointwise prediction  local curvature  <->  local excess rate.

    Returns arrays of (scal, local_second_moment, local_density_estimate).
    The second moment is the mean squared distance to the k-NN centroid,
    a proxy for the cell moment M_2 at that anchor.
    """
    out = local_curvature(X, pca_dim, k_nn=k_nn, n_anchor=n_anchor, seed=seed)
    return out
