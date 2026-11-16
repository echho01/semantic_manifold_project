"""Held-out task projections and paired-encoder quantizer transfer."""
from __future__ import annotations
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans
from rd import split_indices


def linear_probe(X, y, dim, seed=0):
    """Return a rank-dim matrix of shape [dim, ambient_dimension]."""
    X, y = np.asarray(X, float), np.asarray(y)
    _, inverse = np.unique(y, return_inverse=True)
    targets = np.eye(inverse.max() + 1)[inverse]
    targets -= targets.mean(0)
    centred = X - X.mean(0)
    weights = np.linalg.solve(centred.T @ centred + 1e-3 * len(X) * np.eye(X.shape[1]),
                              centred.T @ targets).T
    _, singular, vt = np.linalg.svd(weights, full_matrices=False)
    rank = int(np.sum(singular > max(singular[0], 1e-30) * 1e-10))
    if not 1 <= dim <= rank:
        raise ValueError(f"probe_dim must lie between 1 and fitted task rank {rank}")
    return singular[:dim, None] * vt[:dim]


def task_rd(X, y, probe_dim, sizes, n_fit=100000, seed=0):
    X, y = np.asarray(X, float), np.asarray(y)
    train_idx, test_idx = split_indices(len(X), seed)
    test_idx = test_idx[:n_fit]
    train, test = X[train_idx], X[test_idx]
    weights = linear_probe(train, y[train_idx], probe_dim, seed)
    mean = train.mean(0)
    z_train, z_test = (train - mean) @ weights.T, (test - mean) @ weights.T
    out = {"sizes": [], "D_task": [], "probe_dim": int(weights.shape[0]),
           "train_count": len(train), "test_count": len(test), "skipped_sizes": []}
    for n in sizes:
        if n < 1 or n > len(train):
            out["skipped_sizes"].append(int(n))
            continue
        km = MiniBatchKMeans(n_clusters=n, batch_size=max(1024, n), n_init=3,
                             random_state=seed).fit(z_train)
        d2 = np.sum((z_test - km.cluster_centers_[km.predict(z_test)]) ** 2, axis=1)
        out["sizes"].append(int(n))
        out["D_task"].append(float(d2.mean()))
    return out


def task_exponent(curve):
    sizes, dist = np.asarray(curve["sizes"], float), np.asarray(curve["D_task"], float)
    valid = (sizes > 0) & (dist > 0) & np.isfinite(dist)
    if valid.sum() < 2:
        return None
    slope = np.polyfit(np.log(sizes[valid]), np.log(dist[valid]), 1)[0]
    return float(-2 / slope) if slope < 0 else None


def procrustes(A, B):
    """Equal-dimension orthogonal rotation and translation of B onto A."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    if A.shape != B.shape:
        raise ValueError("orthogonal Procrustes requires paired equal-dimensional matrices")
    ac, bc = A - A.mean(0), B - B.mean(0)
    u, _, vt = np.linalg.svd(bc.T @ ac, full_matrices=False)
    return bc @ (u @ vt) + A.mean(0)


def estimate_epsilon(A, B, k=16, n=20000, seed=0):
    """Squared-distance deviations on the SAME paired A-neighbour edges."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    if len(A) != len(B):
        raise ValueError("encoder samples must be paired")
    idx = np.random.default_rng(seed).choice(len(A), min(n, len(A)), replace=False)
    A, B = A[idx], B[idx]
    k = min(k, len(A) - 1)
    _, neighbours = NearestNeighbors(n_neighbors=k + 1).fit(A).kneighbors(A)
    deviations = []
    for i, ids in enumerate(neighbours):
        ids = ids[ids != i][:k]
        da = np.sum((A[ids] - A[i]) ** 2, axis=1)
        db = np.sum((B[ids] - B[i]) ** 2, axis=1)
        valid = da > 1e-20
        deviations.extend(np.abs(db[valid] / da[valid] - 1).tolist())
    if not deviations:
        raise ValueError("no positive-distance paired edges")
    values = np.asarray(deviations)
    return {"median": float(np.median(values)), "q95": float(np.quantile(values, .95)),
            "max_observed": float(values.max()), "n_pairs": len(values),
            "interpretation": "finite-sample paired-edge diagnostic; not a uniform metric bound"}


def mismatch_rd(A, B, sizes, n_fit=50000, seed=0):
    """Preserve A's encoded index and decode to the corresponding B centroid.

    Both the transferred quantizer and its matched baseline are evaluated
    under B's distortion on the identical held-out pairs.
    """
    A, B = np.asarray(A, float), np.asarray(B, float)
    if len(A) != len(B):
        raise ValueError("encoder samples must be paired")
    train_idx, test_idx = split_indices(len(A), seed)
    test_idx = test_idx[:n_fit]
    at, bt, ae, be = A[train_idx], B[train_idx], A[test_idx], B[test_idx]
    out = {"sizes": [], "D_A_matched": [], "D_B_matched": [], "D_B_transferred": [],
           "transfer_over_matched_B": [], "skipped_sizes": [],
           "train_count": len(train_idx), "test_count": len(test_idx)}
    for n in sizes:
        if n < 1 or n > len(at):
            out["skipped_sizes"].append(int(n))
            continue
        params = dict(n_clusters=n, batch_size=max(1024, n), n_init=3, random_state=seed)
        ka, kb = MiniBatchKMeans(**params).fit(at), MiniBatchKMeans(**params).fit(bt)
        labels = ka.predict(at)
        mapped = np.zeros((n, B.shape[1]))
        mass = np.bincount(labels, minlength=n)
        np.add.at(mapped, labels, bt)
        good = mass > 0
        mapped[good] /= mass[good, None]
        # Rare empty learned cells map through their nearest training pair.
        if not good.all():
            nearest = NearestNeighbors(n_neighbors=1).fit(at).kneighbors(ka.cluster_centers_[~good], return_distance=False)[:, 0]
            mapped[~good] = bt[nearest]
        encoded = ka.predict(ae)
        da = float(np.mean(np.sum((ae - ka.cluster_centers_[encoded]) ** 2, axis=1)))
        db = float(np.mean(np.sum((be - kb.cluster_centers_[kb.predict(be)]) ** 2, axis=1)))
        transferred = float(np.mean(np.sum((be - mapped[encoded]) ** 2, axis=1)))
        out["sizes"].append(int(n))
        out["D_A_matched"].append(da)
        out["D_B_matched"].append(db)
        out["D_B_transferred"].append(transferred)
        out["transfer_over_matched_B"].append(transferred / max(db, 1e-30))
    return out
