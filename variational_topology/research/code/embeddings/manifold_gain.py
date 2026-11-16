"""Manifold gain and semantic outage on real embeddings.

Two experiments:

1. Gain (Theorem 4.1 / Proposition 4.2).
   Add Gaussian noise of total power sigma^2 in ambient dimension D_a to
   held-out embeddings, project back onto the estimated tangent plane at
   each point (local PCA with k neighbours), and measure

       geo_mse   = mean squared tangential displacement (the distortion)
       raw_mse   = mean squared ambient displacement (= sigma^2)
       gain_dB   = 10 log10(raw_mse / geo_mse)

   The prediction is gain_dB -> 10 log10(D_a / d) for sigma small compared
   with the local radius, with a curvature-dependent O(sigma^2) correction
   whose sign is given by Proposition 4.2 on space forms.

2. Outage (Theorem 4.4).
   Classification by nearest class centroid.  The concept separation is
   Delta_0 = 5th percentile of the pairwise distance between different
   classes' centroids; the prediction is that the error curve is flat
   for sigma << min(tau, Delta_0) and rises at sigma ~ min(tau, Delta_0),
   with the reach estimated from the local curvature fit.
"""
from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


def tangent_project(y: np.ndarray, X: np.ndarray, k: int = 32) -> np.ndarray:
    """Project each row of y onto the local PCA tangent plane of X at its
    nearest neighbour."""
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    _, idx = nn.kneighbors(y)
    idx = idx[:, 1:]
    out = np.empty_like(y)
    for i in range(len(y)):
        V = X[idx[i]] - X[idx[i]].mean(0)
        _, _, Vt = np.linalg.svd(V, full_matrices=False)
        # tangent basis: the top-k singular directions of the local cloud
        E = Vt[: V.shape[0] // 2 or 1].T
        p = y[i] @ E
        out[i] = p @ E.T
    return out


def gain_curve(X: np.ndarray, levels: list[float], n_points: int = 20000,
               k: int = 32, seed: int = 0) -> dict:
    """Projection decoder: orthogonal projection of the noisy point onto the
    local PCA plane at its *un-noised* anchor.

    Using the known anchor isolates the effect the theorem describes (the
    noise power that survives projection), rather than the anchor-search
    error, which is a separate and much larger effect.  Operationally this
    models a receiver that knows the codeword and projects the observation
    onto the codeword's tangent plane -- the projection decoder of
    Section 4.
    """
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), min(n_points, len(X)), replace=False)
    Xs = X[idx].astype(np.float32)
    n, D = Xs.shape
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Xs)
    _, nbr = nn.kneighbors(Xs)
    nbr = nbr[:, 1:]
    bases = []
    for i in range(len(Xs)):
        V = Xs[nbr[i]] - Xs[nbr[i]].mean(0)
        _, _, Vt = np.linalg.svd(V, full_matrices=False)
        bases.append(Vt.T)                      # [D, k]
    out = {"sigma": [], "gain_dB": [], "geo_mse": [], "raw_mse": []}
    for s_tot in levels:
        s_ax = s_tot / np.sqrt(D)
        Y = Xs + rng.normal(scale=s_ax, size=Xs.shape).astype(np.float32)
        geo = np.empty(len(Y), np.float32)
        raw = np.empty(len(Y), np.float32)
        for i in range(len(Y)):
            B = bases[i]
            d = Y[i] - Xs[i]
            P = B @ (B.T @ d)
            geo[i] = float((P ** 2).sum())
            raw[i] = float((d ** 2).sum())
        out["sigma"].append(float(s_tot))
        out["geo_mse"].append(float(geo.mean()))
        out["raw_mse"].append(float(raw.mean()))
        out["gain_dB"].append(float(10 * np.log10(raw.mean() / max(geo.mean(), 1e-30))))
    return out


def class_separation(X: np.ndarray, labels: np.ndarray) -> float:
    """5th percentile of the distance between different class centroids."""
    labs = np.unique(labels)
    C = np.stack([X[labels == l].mean(0) for l in labs])
    D2 = ((C[:, None] - C[None]) ** 2).sum(-1)
    iu = np.triu_indices(len(C), 1)
    return float(np.sqrt(np.percentile(D2[iu], 5)))


def outage_curve(X: np.ndarray, labels: np.ndarray, levels: list[float],
                 n_points: int = 20000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), min(n_points, len(X)), replace=False)
    Xs, ys = X[idx].astype(np.float32), labels[idx]
    D = Xs.shape[1]
    labs = np.unique(ys)
    C = np.stack([Xs[ys == l].mean(0) for l in labs])
    out = {"sigma": [], "error": [], "m": []}
    delta0 = class_separation(Xs, ys)
    for s_tot in levels:
        s_ax = s_tot / np.sqrt(D)
        Y = Xs + rng.normal(scale=s_ax, size=Xs.shape).astype(np.float32)
        pred = np.argmin(((Y[:, None, :] - C[None]) ** 2).sum(-1), 1)
        true = np.searchsorted(labs, ys)
        out["sigma"].append(float(s_tot))
        out["error"].append(float((pred != true).mean()))
        out["m"].append(delta0)
    out["delta0"] = delta0
    return out
