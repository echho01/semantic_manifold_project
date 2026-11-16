"""Persistent-homology utilities shared by the Paper 2 scripts.

Conventions (match paper2-topology/main.tex, Section II):

  * Diagrams are in the RADIUS convention.  A Rips simplex enters at half
    its diameter, so H0 deaths are half the single-linkage merge heights
    and Rips/Cech agree in degree 0.  ripser reports edge lengths; every
    wrapper below divides by 2.
  * Bars are half-open [b, d); essential classes have d = inf.
  * Noise is normalized as in Paper 1: N ~ N(0, sigma^2/D_a I), so that
    E||N||^2 = sigma^2 and ||N_i - N_j||^2 ~ 2 sigma^2.

Functions
  rips             Rips diagrams of a point cloud or a distance matrix
  h0_deaths        exact H0 deaths from the minimum spanning tree
  bottleneck       exact bottleneck distance (finite candidate search + matching)
  n_significant    number of bars with persistence > theta
  rank_estimator   rank of H(F_a) -> H(F_b): bars with b_i <= a, d_i > b
  phi_reparam      the high-dimensional reparametrization r -> sqrt(r^2+s^2/2)
  debias           distance debiasing sqrt((d^2 - 2 sigma^2)_+)
  maxmin           greedy farthest-point landmarks
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching, minimum_spanning_tree
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform

RIPS_SANDWICH = 2.0 + np.sqrt(2.0)   # Rips inference uses scales (eta, (2+sqrt2) eta)
RIPS_REACH = 3.0 + 2.0 * np.sqrt(2.0)  # ... valid when RIPS_REACH * eta < tau


# ------------------------------------------------------------------ diagrams
def _distance_matrix(D: np.ndarray) -> np.ndarray:
    D = np.asarray(D, dtype=np.float64)
    if D.ndim != 2 or D.shape[0] != D.shape[1] or not len(D):
        raise ValueError("D must be a nonempty square dissimilarity matrix")
    if not np.isfinite(D).all() or np.any(D < 0):
        raise ValueError("dissimilarities must be finite and nonnegative")
    if not np.allclose(D, D.T, rtol=1e-12, atol=1e-14):
        raise ValueError("dissimilarities must be symmetric")
    if np.any(np.diag(D) != 0):
        raise ValueError("dissimilarity diagonal must be zero")
    return D


def rips(X: np.ndarray | None = None, maxdim: int = 1, D: np.ndarray | None = None,
         thresh: float = np.inf) -> list[np.ndarray]:
    """Rips diagrams in the radius convention, over F_2.

    D may be a dissimilarity matrix: the triangle inequality is unnecessary
    for a flag filtration. With a finite thresh, infinite deaths are
    right-censored intervals, not evidence of globally essential classes.
    """
    from ripser import ripser
    if (X is None) == (D is None):
        raise ValueError("provide exactly one of X and D")
    if int(maxdim) != maxdim or maxdim < 0 or thresh < 0:
        raise ValueError("invalid homology degree or filtration cutoff")
    if D is not None:
        res = ripser(_distance_matrix(D), maxdim=maxdim, distance_matrix=True,
                     thresh=2 * thresh, coeff=2)
    else:
        X = np.asarray(X, dtype=np.float64)
        if X.ndim != 2 or not len(X) or not np.isfinite(X).all():
            raise ValueError("X must be a nonempty finite point cloud")
        res = ripser(X, maxdim=maxdim, thresh=2 * thresh, coeff=2)
    return [np.asarray(d, float).reshape(-1, 2) / 2.0 for d in res["dgms"]]


def h0_deaths(X: np.ndarray | None = None, D: np.ndarray | None = None) -> np.ndarray:
    """All n-1 finite H0 death radii, including zero deaths at duplicates.

    Sparse graph libraries reserve zero entries for missing edges. Replace
    zero-weight off-diagonal edges by the smallest positive float while
    constructing the MST, then restore their true zero weight.
    """
    if (X is None) == (D is None):
        raise ValueError("provide exactly one of X and D")
    if D is None:
        D = squareform(pdist(np.asarray(X, float)))
    D = _distance_matrix(D)
    n = len(D)
    if n == 1:
        return np.empty(0)
    i, j = np.triu_indices(n, 1)
    tiny = np.nextafter(0.0, 1.0)
    weights = D[i, j].copy()
    weights[weights == 0] = tiny
    T = minimum_spanning_tree(csr_matrix((weights, (i, j)), shape=(n, n)))
    deaths = T.data.copy()
    deaths[deaths == tiny] = 0.0
    if len(deaths) != n - 1:
        raise RuntimeError("complete dissimilarity graph has no spanning tree")
    return np.sort(deaths)[::-1] / 2.0


def h0_diagram(deaths: np.ndarray) -> np.ndarray:
    """H0 diagram (with one essential class) from finite deaths."""
    deaths = np.asarray(deaths, float)
    return np.vstack([np.c_[np.zeros_like(deaths), deaths], [[0.0, np.inf]]])


def persistence(dgm: np.ndarray) -> np.ndarray:
    return dgm[:, 1] - dgm[:, 0]


def n_significant(dgm: np.ndarray, theta: float) -> int:
    """Bars with persistence > theta (essential bars always count)."""
    if len(dgm) == 0:
        return 0
    return int(np.sum(persistence(dgm) > theta))


def rank_estimator(dgm: np.ndarray, a: float, b: float) -> int:
    """rank H(F_a) -> H(F_b) = #{bars [b_i, d_i) with b_i <= a < b < d_i}."""
    if len(dgm) == 0:
        return 0
    return int(np.sum((dgm[:, 0] <= a) & (dgm[:, 1] > b)))


# ---------------------------------------------------------------- bottleneck
def _split(dgm: np.ndarray):
    dgm = np.asarray(dgm, float).reshape(-1, 2)
    if np.any(~np.isfinite(dgm[:, 0])) or np.any(np.isnan(dgm[:, 1])) or np.any(dgm[:, 1] < dgm[:, 0]):
        raise ValueError("diagram intervals must have finite births and birth <= death")
    fin = np.isfinite(dgm[:, 1])
    F = dgm[fin]
    F = F[F[:, 1] > F[:, 0]]            # drop zero-length bars
    return F, np.sort(dgm[~fin, 0])


def _covers(P: np.ndarray, Q: np.ndarray, heavyP: np.ndarray, t: float) -> bool:
    """Is there a matching in the L_inf-threshold-t graph P x Q that covers
    every heavy point of P?"""
    hp = np.flatnonzero(heavyP)
    if len(hp) == 0:
        return True
    if len(Q) < len(hp):
        return False
    tree = cKDTree(Q)
    rows, cols = [], []
    for r, i in enumerate(hp):
        nb = tree.query_ball_point(P[i], np.nextafter(t, np.inf), p=np.inf)
        nb = [j for j in nb if np.max(np.abs(P[i] - Q[j])) <= t]
        if not nb:
            return False
        rows.extend([r] * len(nb))
        cols.extend(nb)
    G = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                   shape=(len(hp), len(Q)))
    m = maximum_bipartite_matching(G, perm_type="column")
    return bool(np.all(m >= 0))


def bottleneck(A: np.ndarray, B: np.ndarray, tol: float | None = None) -> float:
    """Bottleneck distance by exact search of floating-point candidate costs.

    Each feasibility test asks for a matching covering all points unable
    to reach the diagonal. The two one-sided saturation conditions are
    equivalent to simultaneous saturation (Dulmage--Mendelsohn). ``tol`` is
    retained for compatibility; no bisection tolerance affects the result.
    """
    FA, EA = _split(A)
    FB, EB = _split(B)
    if len(EA) != len(EB):
        return np.inf
    ess = float(np.max(np.abs(EA - EB))) if len(EA) else 0.0
    pa = (FA[:, 1] - FA[:, 0]) / 2 if len(FA) else np.empty(0)
    pb = (FB[:, 1] - FB[:, 0]) / 2 if len(FB) else np.empty(0)
    if not len(FA) or not len(FB):
        return max(ess, float(pa.max(initial=0)), float(pb.max(initial=0)))
    cross = np.max(np.abs(FA[:, None, :] - FB[None, :, :]), axis=2)
    upper = max(float(pa.max(initial=0)), float(pb.max(initial=0)))
    candidates = np.unique(np.r_[0.0, pa, pb, cross[cross <= upper]])

    def saturates(costs, heavy, t):
        selected = costs[heavy]
        if len(selected) == 0:
            return True
        if len(selected) > costs.shape[1]:
            return False
        graph = csr_matrix(selected <= t)
        matched = maximum_bipartite_matching(graph, perm_type="column")
        return bool(np.all(matched >= 0))

    def feasible(t):
        return saturates(cross, pa > t, t) and saturates(cross.T, pb > t, t)

    lo, hi = 0, len(candidates) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if feasible(candidates[mid]):
            hi = mid
        else:
            lo = mid + 1
    return max(float(candidates[lo]), ess)


def bottleneck_by_degree(dA: list, dB: list) -> list[float]:
    if len(dA) != len(dB):
        raise ValueError("compare the same homology degrees")
    return [bottleneck(a, b) for a, b in zip(dA, dB)]


# --------------------------------------------------- high-dimensional noise
def phi(r: np.ndarray, sigma: float) -> np.ndarray:
    """Radius-convention reparametrization induced by isotropic noise of
    total power sigma^2:  edge l -> sqrt(l^2 + 2 sigma^2)."""
    return np.sqrt(np.asarray(r, float) ** 2 + sigma ** 2 / 2.0)


def noise_reference(D: np.ndarray, sigma: float) -> np.ndarray:
    """Deterministic noise-reference dissimilarities, including duplicate pairs.

    Off-diagonal entries are sqrt(D_ij**2 + 2*sigma**2); vertices enter at
    zero. This is the safe reference when a clean persistence diagram has
    discarded zero-length H0 bars.
    """
    D = _distance_matrix(D)
    if sigma < 0 or not np.isfinite(sigma):
        raise ValueError("sigma must be finite and nonnegative")
    out = np.sqrt(D * D + 2 * sigma * sigma)
    np.fill_diagonal(out, 0)
    return out


def phi_reparam(dgm: np.ndarray, sigma: float, degree: int) -> np.ndarray:
    """Transform bar endpoints; H0 input must include all zero deaths.

    Use ``noise_reference`` when diagrams came from a library which omits
    diagonal H0 bars. Vertices always retain birth zero.
    """
    out = np.array(dgm, float, copy=True).reshape(-1, 2)
    if degree > 0:
        out[:, 0] = phi(out[:, 0], sigma)
    fin = np.isfinite(out[:, 1])
    out[fin, 1] = phi(out[fin, 1], sigma)
    return out


def sq_dists(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, np.float64)
    g = (X * X).sum(1)
    D2 = g[:, None] + g[None, :] - 2.0 * (X @ X.T)
    np.maximum(D2, 0.0, out=D2)
    np.fill_diagonal(D2, 0.0)
    return D2


def debias(X: np.ndarray, sigma: float) -> np.ndarray:
    """Debiased distance matrix sqrt((||x_i-x_j||^2 - 2 sigma^2)_+)."""
    D2 = sq_dists(X) - 2.0 * sigma ** 2
    np.maximum(D2, 0.0, out=D2)
    np.fill_diagonal(D2, 0.0)
    return np.sqrt(D2)


def zeta_raw(sigma: float, n: int, Da: int, delta: float = 0.05) -> float:
    """Theorem bound on d_B(noisy, phi_* clean), radius convention."""
    a = np.log(2 * n * n / delta) / Da
    return sigma * ((2 + np.sqrt(2)) * np.sqrt(a) + np.sqrt(2) * a)


def zeta_debiased(sigma: float, n: int, Da: int, delta: float = 0.05) -> float:
    """Theorem bound on d_B(debiased, clean), radius convention."""
    a = np.sqrt(np.log(2 * n * n / delta) / Da)
    alpha, beta = 4 * sigma * a, 4 * sigma ** 2 * (a + a * a)
    return (alpha + np.sqrt(alpha ** 2 + 4 * beta)) / 4.0


# ------------------------------------------------------- Cech (alpha) in R^2
def alpha2d(X: np.ndarray, backend: str = "auto") -> list[np.ndarray]:
    """Cech-equivalent (alpha complex) diagrams of a planar point set, radius
    convention, degrees 0 and 1.  Uses gudhi when available; otherwise a
    Delaunay construction with the Gabriel rule and Z/2 column reduction.
    """
    X = np.unique(np.asarray(X, float), axis=0)
    if X.ndim != 2 or X.shape[1] != 2 or not len(X) or not np.isfinite(X).all():
        raise ValueError("alpha2d requires a nonempty finite planar cloud")
    if backend not in {"auto", "gudhi", "native"}:
        raise ValueError("unknown alpha backend")
    if len(X) < 3 or np.linalg.matrix_rank(X - X.mean(0)) < 2:
        return [h0_diagram(h0_deaths(X=X)), np.empty((0, 2))]
    try:
        if backend == "native":
            raise ImportError
        import gudhi
        st = gudhi.AlphaComplex(points=X).create_simplex_tree()
        st.compute_persistence()
        out = []
        for k in (0, 1):
            I = st.persistence_intervals_in_dimension(k)
            I = np.sqrt(np.maximum(np.asarray(I, float).reshape(-1, 2), 0.0))
            out.append(I)
        return out
    except ImportError:
        if backend == "gudhi":
            raise
    from scipy.spatial import Delaunay
    tri = Delaunay(X).simplices
    n = len(X)
    # triangle circumradii
    A, B, C = X[tri[:, 0]], X[tri[:, 1]], X[tri[:, 2]]
    a = np.linalg.norm(B - C, axis=1)
    b = np.linalg.norm(A - C, axis=1)
    c = np.linalg.norm(A - B, axis=1)
    area2 = np.abs((B[:, 0] - A[:, 0]) * (C[:, 1] - A[:, 1])
                   - (B[:, 1] - A[:, 1]) * (C[:, 0] - A[:, 0]))
    rho = a * b * c / (2.0 * np.maximum(area2, 1e-300))
    # edges with incident triangles and opposite vertices
    edge_id, edge_inc = {}, []
    tri_edges = np.empty((len(tri), 3), np.int64)
    for t, (i, j, k) in enumerate(tri):
        for s, (u, v, w) in enumerate(((i, j, k), (j, k, i), (i, k, j))):
            key = (u, v) if u < v else (v, u)
            if key not in edge_id:
                edge_id[key] = len(edge_inc)
                edge_inc.append([])
            e = edge_id[key]
            edge_inc[e].append((t, w))
            tri_edges[t, s] = e
    edges = np.array(list(edge_id.keys()), np.int64)
    fe = np.empty(len(edges))
    for e, (u, v) in enumerate(edges):
        m = 0.5 * (X[u] + X[v])
        r = 0.5 * np.linalg.norm(X[u] - X[v])
        gabriel = all(np.linalg.norm(X[w] - m) >= r for _, w in edge_inc[e])
        fe[e] = r if gabriel else min(rho[t] for t, _ in edge_inc[e])
    ft = np.maximum(rho, fe[tri_edges].max(1))      # monotone by construction
    # H0 by Kruskal
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    order_e = np.argsort(fe, kind="stable")
    rank_e = np.empty(len(edges), np.int64)
    rank_e[order_e] = np.arange(len(edges))
    h0, positive = [], np.zeros(len(edges), bool)
    for e in order_e:
        ru, rv = find(edges[e, 0]), find(edges[e, 1])
        if ru == rv:
            positive[e] = True
        else:
            parent[ru] = rv
            h0.append((0.0, fe[e]))
    h0.append((0.0, np.inf))
    # H1 by column reduction of triangle boundaries (pivot = youngest edge)
    pivot_of = {}
    h1 = []
    for t in np.argsort(ft, kind="stable"):
        col = set(int(rank_e[e]) for e in tri_edges[t] if positive[e])
        while col:
            low = max(col)
            if low in pivot_of:
                col ^= pivot_of[low]
            else:
                pivot_of[low] = col
                h1.append((fe[order_e[low]], ft[t]))
                break
    h1 = np.asarray(h1).reshape(-1, 2)
    h1 = h1[h1[:, 1] > h1[:, 0]]
    return [np.asarray(h0), h1]


# ---------------------------------------------------------------- sampling
def maxmin(X: np.ndarray, m: int, seed: int = 0) -> tuple[np.ndarray, float]:
    """Greedy farthest-point landmarks; returns (indices, covering radius)."""
    rng = np.random.default_rng(seed)
    X = np.asarray(X, float)
    n = len(X)
    if not n or m < 1 or int(m) != m:
        raise ValueError("landmark count must be a positive integer")
    m = min(m, n)
    chosen = np.zeros(n, bool)
    idx = np.empty(m, np.int64)
    idx[0] = rng.integers(n)
    chosen[idx[0]] = True
    d = np.linalg.norm(X - X[idx[0]], axis=1)
    for i in range(1, m):
        idx[i] = int(np.argmax(np.where(chosen, -np.inf, d)))
        chosen[idx[i]] = True
        d = np.minimum(d, np.linalg.norm(X - X[idx[i]], axis=1))
    return idx, float(d.max())


def random_frame(dim_in: int, Da: int, seed: int = 0) -> np.ndarray:
    """Random orthonormal dim_in-frame in R^Da (isometric embedding)."""
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(Da, dim_in)))
    return Q.T                         # [dim_in, Da]
