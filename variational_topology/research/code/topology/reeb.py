"""Reeb-graph tools for task functionals (Paper 2, Section VI).

  level_components   connected components of the slabs f^{-1}[t_j, t_{j+1})
                     by single linkage at scale `link`; returns a slab index
                     and a component label per sample
  branch_entropy     plug-in estimate of H(E | f(S)) in bits: the entropy of
                     the level-set component given the (binned) task value
  mapper_graph       Mapper graph (overlapping slabs), nodes and edges, with
                     its first Betti number  #E - #V + #components
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def _components(X: np.ndarray, link: float) -> np.ndarray:
    if len(X) == 0:
        return np.zeros(0, np.int64)
    pairs = cKDTree(X).query_pairs(link, output_type="ndarray")
    from scipy.sparse import coo_matrix
    n = len(X)
    A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(A, directed=False)[1]


def level_components(X: np.ndarray, f: np.ndarray, n_bins: int, link: float,
                     min_size: int = 1):
    edges = np.quantile(f, np.linspace(0, 1, n_bins + 1))
    slab = np.clip(np.searchsorted(edges, f, side="right") - 1, 0, n_bins - 1)
    comp = np.full(len(f), -1, np.int64)
    for j in range(n_bins):
        m = np.flatnonzero(slab == j)
        lab = _components(X[m], link)
        # merge tiny components (noise) into their nearest large one
        sizes = np.bincount(lab) if len(lab) else np.zeros(0, int)
        big = np.flatnonzero(sizes >= min_size)
        if len(big) and len(big) < len(sizes):
            cent = {c: X[m][lab == c] for c in big}
            for c in np.flatnonzero(sizes < min_size):
                pts = X[m][lab == c]
                dist = [np.min(np.linalg.norm(cent[b][:, None] - pts[None], axis=2)) for b in big]
                lab[lab == c] = big[int(np.argmin(dist))]
        comp[m] = lab
    return slab, comp


def branch_entropy(X: np.ndarray, f: np.ndarray, n_bins: int, link: float,
                   min_size: int = 1) -> float:
    """H(E | f) ~ sum_j P(slab j) H(component | slab j), bits."""
    slab, comp = level_components(X, f, n_bins, link, min_size)
    H = 0.0
    for j in np.unique(slab):
        c = comp[slab == j]
        p = np.bincount(c - c.min()) / len(c)
        p = p[p > 0]
        H += (len(c) / len(f)) * float(-(p * np.log2(p)).sum())
    return H


def mapper_graph(X: np.ndarray, f: np.ndarray, n_int: int, overlap: float,
                 link: float, min_size: int = 1) -> dict:
    if n_int < 2 or not 0 < overlap < 0.5:
        raise ValueError("use at least two intervals with 0 < overlap < 0.5 so the nerve is a graph")
    lo, hi = float(f.min()), float(f.max())
    if hi <= lo:
        lab = _components(X, link)
        count = int(sum(np.sum(lab == c) >= min_size for c in np.unique(lab)))
        return {"nodes": count, "edges": 0, "beta0": count, "beta1": 0}
    w = (hi - lo) / (n_int - (n_int - 1) * overlap)
    nodes = []
    for j in range(n_int):
        a = lo + j * w * (1 - overlap)
        m = np.flatnonzero((f >= a) & (f <= a + w))
        lab = _components(X[m], link)
        for c in np.unique(lab):
            mem = m[lab == c]
            if len(mem) >= min_size:
                nodes.append(set(mem.tolist()))
    E = set()
    for i in range(len(nodes)):
        for k in range(i + 1, len(nodes)):
            if nodes[i] & nodes[k]:
                E.add((i, k))
    V = len(nodes)
    if V == 0:
        return {"nodes": 0, "edges": 0, "beta0": 0, "beta1": 0}
    from scipy.sparse import coo_matrix
    if E:
        e = np.array(sorted(E))
        A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(V, V))
    else:
        A = coo_matrix((V, V))
    b0 = connected_components(A, directed=False)[0]
    return {"nodes": V, "edges": len(E), "beta0": int(b0),
            "beta1": int(len(E) - V + b0)}
