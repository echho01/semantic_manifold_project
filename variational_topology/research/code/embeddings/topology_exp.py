"""Paper 2 (topology) measurements on real embeddings.

Reads the arrays.npz written by the Paper 1 extract stage
(runs/<model>__<dataset>/arrays.npz) and implements the protocol of
paper2-topology/EXPERIMENTS.md.  Every function returns plain dicts.

  common_items     fixed item subset shared by all models (same seed)
  normalize        scale so the median pairwise distance of the subset is 1
  diagrams         Rips H0/H1 (radius convention) of a cloud
  floor_matrix     d_B between models on the same items        (Thm 19)
  dictionary_curve d_topo(K) with dictionaries chosen in model T and
                   interpreted in model R, plus eps_R(K)          (Thm 19)
  reparam_check    d_B(noisy, phi_* clean) <= zeta_raw on real data (Thm 11)
  cliff_curve      raw / debiased / KB-aware receivers vs sigma  (Thm 14, Cor 16)
  reeb_overhead    binned component entropy and Mapper graph diagnostics
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "topology"))
import reeb  # noqa: E402
import tda  # noqa: E402


# ------------------------------------------------------------- utilities
def common_items(n_total: int, n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_total, size=min(n, n_total), replace=False))


def normalization(X: np.ndarray, seed: int = 0, n_ref: int = 1000):
    """Fit centering and scale on a reference (training) split only."""
    X = np.asarray(X, np.float64)
    if X.ndim != 2 or len(X) < 2 or not np.isfinite(X).all():
        raise ValueError("normalization needs at least two finite feature rows")
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), size=min(n_ref, len(X)), replace=False)
    d = np.sqrt(tda.sq_dists(X[idx]))
    med = float(np.median(d[np.triu_indices(len(idx), 1)]))
    if med <= 0:
        raise ValueError("median reference distance is zero")
    return X.mean(0), med


def normalize(X: np.ndarray, seed: int = 0, n_ref: int = 2000) -> np.ndarray:
    """Center and scale so that the median pairwise distance is 1.  Makes
    diagrams of different models comparable; d_B is then in units of the
    model's own typical distance."""
    mean, scale = normalization(X, seed, n_ref)
    return (np.asarray(X, np.float64) - mean) / scale


def diagrams(X: np.ndarray, maxdim: int = 1, D: np.ndarray | None = None) -> list:
    return tda.rips(X=None if D is not None else X, D=D, maxdim=maxdim)


def _summary(dgms: list) -> dict:
    out = {}
    for k, d in enumerate(dgms):
        p = tda.persistence(d)
        fin = np.isfinite(p)
        out[f"H{k}"] = {"n_bars": int(len(d)), "n_essential": int((~fin).sum()),
                        "top_persistence": np.sort(p[fin])[::-1][:10].round(5).tolist()}
    return out


def _dB(A: list, B: list) -> float:
    return float(max(tda.bottleneck_by_degree(A, B)))


# ------------------------------------------------------------- floor
def floor_matrix(clouds: dict, maxdim: int = 1) -> dict:
    """clouds: name -> normalized [n, D] array on the SAME items."""
    names = sorted(clouds)
    dg = {m: diagrams(clouds[m], maxdim) for m in names}
    M = {a: {b: (0.0 if a == b else _dB(dg[a], dg[b])) for b in names} for a in names}
    return {"models": names, "floor": M, "diagrams": {m: _summary(dg[m]) for m in names}}


# ------------------------------------------------------------- dictionary
def dictionary_curve(XT: np.ndarray, XR: np.ndarray, sizes: list, maxdim: int = 1,
                     seed: int = 0) -> dict:
    """Dictionary chosen by farthest-point sampling in the transmitter's
    geometry XT; the receiver's structure is the Rips diagram of the same
    items in its own geometry XR.  Reference: Dgm(XT) on all items."""
    if len(XT) != len(XR) or not sizes or min(sizes) < 1 or max(sizes) > len(XT):
        raise ValueError("paired clouds and dictionary sizes must be valid")
    gT, gR = diagrams(XT, maxdim), diagrams(XR, maxdim)
    floor = _dB(gT, gR)
    order, _ = tda.maxmin(XT, max(sizes), seed=seed)
    out = {"floor": floor, "K": [], "d_topo": [], "eps_R": [], "within_bound": []}
    for K in sizes:
        idx = order[:K]
        d = _dB(diagrams(XR[idx], maxdim), gT)
        from scipy.spatial import cKDTree
        epsR = float(cKDTree(XR[idx]).query(XR, k=1)[0].max())
        out["K"].append(int(K))
        out["d_topo"].append(d)
        out["eps_R"].append(epsR)
        out["within_bound"].append(bool(abs(d - floor) <= epsR + 1e-9))
    return out


# ------------------------------------------------------------- reparam
def reparam_check(X: np.ndarray, sigmas: list, maxdim: int = 1, delta: float = 0.05,
                  seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    n, Da = X.shape
    dX = diagrams(X, maxdim)
    clean_distances = np.sqrt(tda.sq_dists(X))
    rows = []
    for s in sigmas:
        Y = X + rng.normal(scale=s / np.sqrt(Da), size=X.shape)
        dY = diagrams(Y, maxdim)
        dP = diagrams(None, maxdim, D=tda.noise_reference(clean_distances, s))
        dD = diagrams(None, maxdim, D=tda.debias(Y, s))
        rows.append({"sigma": float(s),
                     "naive": _dB(dY, dX),
                     "raw_vs_phi": _dB(dY, dP), "zeta_raw": tda.zeta_raw(s, n, Da, delta),
                     "debiased": _dB(dD, dX), "zeta_deb": tda.zeta_debiased(s, n, Da, delta)})
    return {"n": int(n), "D_a": int(Da), "rows": rows,
            "all_within": bool(all(r["raw_vs_phi"] <= r["zeta_raw"] and
                                   r["debiased"] <= r["zeta_deb"] for r in rows))}


# ------------------------------------------------------------- cliffs
def cliff_curve(X: np.ndarray, labels: np.ndarray, sigmas: list, per_class: int = 60,
                theta: float | None = None, seed: int = 0,
                X_reference: np.ndarray | None = None,
                labels_reference: np.ndarray | None = None) -> dict:
    """Held-out receiver comparison; class means and threshold use references.

    The class-count threshold is calibrated on a disjoint reference cloud.
    Correct count and classification accuracy are separate statistics.
    """
    if X_reference is None or labels_reference is None:
        raise ValueError("provide a disjoint reference split for class means and threshold")
    rng = np.random.default_rng(seed)
    classes = np.unique(labels_reference)
    if not np.isin(labels, classes).all():
        raise ValueError("held-out labels include an unseen reference class")

    def balanced_ids(y):
        return np.concatenate([rng.choice(np.flatnonzero(y == c),
            size=min(per_class, int(np.sum(y == c))), replace=False) for c in classes])

    idx = balanced_ids(labels)
    Xc, yc = X[idx], labels[idx]
    ref = X_reference[balanced_ids(labels_reference)]
    n, Da = Xc.shape
    means = np.stack([X_reference[labels_reference == c].mean(0) for c in classes])
    if theta is None:
        deaths = tda.h0_deaths(X=ref)
        k = len(classes) - 1
        if k == 0:
            theta = 1.05 * float(deaths[0])
        elif len(deaths) <= k:
            raise ValueError("reference cloud needs more than one sample per class")
        else:
            theta = float((deaths[k - 1] + deaths[k]) / 2)
    clean = tda.h0_deaths(X=Xc)
    out = {"theta": float(theta), "n_classes": int(len(classes)),
           "threshold_source": "disjoint_reference_cloud", "n_evaluation": n,
           "clean_count": int((clean > theta).sum()) + 1, "rows": []}
    noise = rng.normal(size=Xc.shape) / np.sqrt(Da)
    for s in sigmas:
        Y = Xc + s * noise
        raw = tda.h0_deaths(D=np.sqrt(tda.sq_dists(Y)))
        deb = tda.h0_deaths(D=tda.debias(Y, s))
        d2 = (means ** 2).sum(1)[None, :] - 2.0 * (Y @ means.T)
        acc = float(np.mean(classes[np.argmin(d2, 1)] == yc))
        out["rows"].append({"sigma": float(s),
                            "raw_count": int((raw > theta).sum()) + 1,
                            "debiased_count": int((deb > theta).sum()) + 1,
                            "kb_aware_acc": acc})
    out["pred_spurious_sigma"] = float(np.sqrt(2) * theta)
    return out


# ------------------------------------------------------------- task
def probe_score(X: np.ndarray, y: np.ndarray | None = None, ridge: float = 1e-3,
                X_reference: np.ndarray | None = None,
                labels_reference: np.ndarray | None = None) -> np.ndarray:
    """Evaluate the first ridge-probe score fitted on a disjoint reference split."""
    if X_reference is None or labels_reference is None:
        raise ValueError("a disjoint reference split is required to fit the probe")
    mean = X_reference.mean(0)
    Xc = X_reference - mean
    classes, yi = np.unique(labels_reference, return_inverse=True)
    Y = np.eye(len(classes))[yi]
    W = np.linalg.solve(Xc.T @ Xc + ridge * len(Xc) * np.eye(X.shape[1]), Xc.T @ Y)
    _, _, Vt = np.linalg.svd(W.T, full_matrices=False)
    return (X - mean) @ Vt[0]


def reeb_overhead(X: np.ndarray, f: np.ndarray, n_bins: int = 30, link: float | None = None,
                  n_int: int = 12, overlap: float = 0.3) -> dict:
    if link is None:
        from scipy.spatial import cKDTree
        d, _ = cKDTree(X).query(X, k=11)
        link = float(np.quantile(d[:, -1], 0.9))
    H = reeb.branch_entropy(X, f, n_bins=n_bins, link=link)
    mp = reeb.mapper_graph(X, f, n_int=n_int, overlap=overlap, link=link)
    slab, comp = reeb.level_components(X, f, n_bins, link)
    bmax = int(max(len(np.unique(comp[slab == j])) for j in np.unique(slab)))
    return {"binned_component_entropy_bits": float(H),
            "H_E_given_f_bits": float(H),
            "interpretation": "finite-sample slab-component diagnostic, not a certified Reeb entropy", "log2_max_components": float(np.log2(bmax)),
            "mapper": mp, "link": link}
