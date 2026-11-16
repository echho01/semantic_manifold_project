"""Task-oriented rate-distortion and metric mismatch.

(a) Task (Theorem 5.1).
    A linear probe  f(x) = W (x - mean)  with W of shape [d_f, D] has a
    pullback semi-metric of rank d_f; the task rate-distortion function is
    that of the quotient by the probe fibers.  Operationally: quantize the
    *projected* representation W x and measure the distortion in the same
    projected space.  The prediction is that the exponent is d_f rather
    than d, i.e. the rate at fixed task distortion falls by
    (d - d_f)/2 * log2(1/D) bits.

(b) Mismatch (Theorem 5.2).
    Two encoders A, B on the same items.  After optional Procrustes
    alignment, epsilon_g is estimated from k-NN distance ratios,
        eps_hat = median over pairs of |d_B^2 / d_A^2 - 1|,
    and the quantizer trained under metric A is decoded under metric B.
    The prediction is a multiplicative excess distortion,
        D_true <= (1 + C eps_hat) D_A,
    i.e. the relative excess is flat in D (no floor), in contrast to a
    topological mismatch.
"""
from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans


def linear_probe(X: np.ndarray, y: np.ndarray, dim: int, seed: int = 0) -> np.ndarray:
    """Ridge-regularized least squares probe, reduced to `dim` output
    dimensions by the top singular directions of the weight matrix."""
    mu = X.mean(0)
    Xc = X - mu
    Y = np.eye(int(y.max()) + 1)[y].astype(np.float32)
    W = np.linalg.solve(Xc.T @ Xc + 1e-3 * len(Xc) * np.eye(X.shape[1]),
                        Xc.T @ Y).T                      # [n_classes, D]
    # collapse to `dim` linear combinations carrying the most task variance
    U, S, Vt = np.linalg.svd(W, full_matrices=False)
    return (U[:, :dim] * S[:dim]) @ Vt[:dim]             # [dim, D]


def task_rd(X: np.ndarray, y: np.ndarray, probe_dim: int, sizes: list[int],
            n_fit: int = 100000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    W = linear_probe(X, y, probe_dim, seed)
    Z = (X - X.mean(0)) @ W.T                            # [n, d_f]
    idx = rng.choice(len(Z), min(n_fit, len(Z)), replace=False)
    Zf = Z[idx]
    out = {"sizes": [], "D_task": []}
    for N in sizes:
        km = MiniBatchKMeans(n_clusters=min(N, len(Z)), n_init=3,
                             random_state=seed).fit(Z)
        d2 = ((Zf - km.cluster_centers_[km.predict(Zf)]) ** 2).sum(1)
        out["sizes"].append(int(N))
        out["D_task"].append(float(d2.mean()))
    return out


def task_exponent(curve: dict) -> float:
    """Slope of log D vs log N: the effective dimension is -2/slope."""
    N = np.log(np.asarray(curve["sizes"], float))
    D = np.log(np.asarray(curve["D_task"], float))
    A = np.vstack([N, np.ones_like(N)]).T
    slope = np.linalg.lstsq(A, D, rcond=None)[0][0]
    return float(-2.0 / slope)


# --------------------------------------------------------------- mismatch
def procrustes(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Rotate/scale B onto A (orthogonal Procrustes, no translation)."""
    Ac = A - A.mean(0); Bc = B - B.mean(0)
    U, S, Vt = np.linalg.svd(Bc.T @ Ac, full_matrices=False)
    R = U @ Vt
    return (Bc @ R) + A.mean(0)


def estimate_epsilon(A: np.ndarray, B: np.ndarray, k: int = 16,
                     n: int = 20000, seed: int = 0) -> float:
    """Median relative deviation of squared k-NN distances."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(A), min(n, len(A)), replace=False)
    A, B = A[idx], B[idx]
    dA, _ = NearestNeighbors(n_neighbors=k).fit(A).kneighbors(A)
    dB, _ = NearestNeighbors(n_neighbors=k).fit(B).kneighbors(B)
    rA = (dA ** 2)[:, -1] / np.maximum((dA ** 2)[:, 0], 1e-30)
    rB = (dB ** 2)[:, -1] / np.maximum((dB ** 2)[:, 0], 1e-30)
    return float(np.median(np.abs(rB / rA - 1.0)))


def mismatch_rd(A: np.ndarray, B: np.ndarray, sizes: list[int],
                n_fit: int = 50000, seed: int = 0) -> dict:
    """Train under metric A, decode under metric B; report the ratio."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(A), min(n_fit, len(A)), replace=False)
    Af, Bf = A[idx].astype(np.float32), B[idx].astype(np.float32)
    out = {"sizes": [], "D_matched": [], "D_mismatched": [], "ratio": []}
    for N in sizes:
        km = MiniBatchKMeans(n_clusters=min(N, len(Af)), n_init=3,
                             random_state=seed).fit(Af)
        labA = km.predict(Af)
        Dm = float(((Af - km.cluster_centers_[labA]) ** 2).sum(1).mean())
        # decode each A-point to its codeword, map to B, take B-nearest
        Bmap = Bf[labA]                       # B coords of the chosen codeword's members
        nn = NearestNeighbors(n_neighbors=1).fit(Bf[labA])   # degenerate; see note
        # honest version: use the B-space centroid of each A-cluster
        centsB = np.stack([Bf[labA == c].mean(0) if (labA == c).any() else Bf[0]
                           for c in range(len(km.cluster_centers_))])
        labB = np.argmin(((Bf[:, None, :] - centsB[None]) ** 2).sum(-1), 1) \
            if N <= 512 else None
        if labB is None:                       # chunked nearest-centroid
            labB = np.empty(len(Bf), np.int64)
            for i in range(0, len(Bf), 4096):
                chunk = Bf[i:i + 4096]
                labB[i:i + 4096] = np.argmin(
                    ((chunk[:, None, :] - centsB[None]) ** 2).sum(-1), 1)
        Dx = float(((Bf - centsB[labB]) ** 2).sum(1).mean())
        out["sizes"].append(int(N))
        out["D_matched"].append(Dm)
        out["D_mismatched"].append(Dx)
        out["ratio"].append(Dx / max(Dm, 1e-30))
    return out
