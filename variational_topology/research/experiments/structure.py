"""Exact finite topology for the declared grid and batch filtrations. Pairing is discrete;
gathering active critical values keeps a piecewise differentiable loss in PyTorch."""

from functools import lru_cache
import numpy as np
import torch
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree, breadth_first_order
from scipy.optimize import linear_sum_assignment
import gudhi


@lru_cache(maxsize=16)
def grid(h, w):
    v = np.arange(h * w).reshape(h, w)
    return np.r_[v[:, :-1].ravel(), v[:-1].ravel()], np.r_[
        v[:, 1:].ravel(), v[1:].ravel()
    ]


def query_pairs(img, n, rng, policy="hybrid"):
    """Fixed given the reference and RNG. Never resample pairs from a reconstruction."""
    f = img.detach().cpu().numpy().reshape(-1)
    pairs = rng.integers(len(f), size=(n, 2))
    fg = np.flatnonzero(f >= 0.5)
    if policy == "hybrid" and len(fg) >= 2:
        pairs[: (n + 1) // 2] = rng.choice(fg, size=((n + 1) // 2, 2))
    return pairs


def merge_values(f, pairs):
    h, w = f.shape
    val = f.detach().cpu().contiguous().numpy().ravel()
    i, j = grid(h, w)
    # Constant positive shift avoids scipy dropping zero-weight edges and leaves MST order unchanged.
    weight = np.maximum(val[i], val[j]).astype(np.float64) + 1.0
    tree = minimum_spanning_tree(
        coo_matrix((weight, (i, j)), shape=(h * w, h * w)).tocsr()
    )
    order, parent = breadth_first_order(tree + tree.T, 0, directed=False)
    depth = np.zeros(h * w, dtype=int)
    for v in order[1:]:
        depth[v] = depth[parent[v]] + 1
    critical = []
    for a, b in pairs:
        u, v = int(a), int(b)
        c = u
        while u != v:
            if depth[u] >= depth[v]:
                if val[u] > val[c]:
                    c = u
                u = parent[u]
            else:
                if val[v] > val[c]:
                    c = v
                v = parent[v]
        if val[u] > val[c]:
            c = u
        critical.append(c)
    return f.reshape(-1)[torch.as_tensor(critical, device=f.device)]


def connection(x, y, n, rng, policy="hybrid"):
    terms = []
    for a, b in zip(x[:, 0], y[:, 0]):
        pairs = query_pairs(a, n, rng, policy)
        terms.append(
            (merge_values(1 - a, pairs) - merge_values(1 - b, pairs)).abs().mean()
        )
    return torch.stack(terms).mean()


def cubical(f):
    # GUDHI takes Fortran-order flattened top-dimensional cells.
    values = f.T.contiguous().reshape(-1)
    cc = gudhi.CubicalComplex(
        dimensions=list(f.shape),
        top_dimensional_cells=values.detach().cpu().double().numpy(),
    )
    cc.persistence(homology_coeff_field=2, min_persistence=0)
    pairs, _ = cc.cofaces_of_persistence_pairs()
    return [
        values[
            torch.as_tensor(pairs[k], dtype=torch.long, device=f.device).reshape(-1, 2)
        ]
        if k < len(pairs)
        else values.new_empty((0, 2))
        for k in (0, 1)
    ]


def rips(d):
    st = gudhi.RipsComplex(
        distance_matrix=d.detach().cpu().double().tolist()
    ).create_simplex_tree(max_dimension=2)
    st.persistence(homology_coeff_field=2, min_persistence=0, persistence_dim_max=True)
    dims = [[], []]

    def critical(simplex):
        if len(simplex) == 1:
            return d.new_zeros(())
        pairs = [(a, b) for l, a in enumerate(simplex) for b in simplex[l + 1 :]]
        vals = torch.stack([d[a, b] for a, b in pairs])
        return vals.max()

    for birth, death in st.persistence_pairs():
        k = len(birth) - 1
        if k in (0, 1) and death:
            dims[k].append(torch.stack((critical(birth), critical(death))))
    return [torch.stack(z) if z else d.new_empty((0, 2)) for z in dims]


def w2(a, b):
    """Wasserstein order 2, infinity ground norm, unlimited matches to the diagonal."""
    n, m = len(a), len(b)
    if not n + m:
        return (a.sum() + b.sum()) * 0
    cost = a.new_full((n + m, m + n), float("inf"))
    if n and m:
        cost[:n, :m] = (a[:, None] - b).abs().amax(2).square()
    i = torch.arange(n, device=a.device)
    j = torch.arange(m, device=a.device)
    cost[i, m + i] = ((a[:, 1] - a[:, 0]) / 2).square()
    cost[n + j, j] = ((b[:, 1] - b[:, 0]) / 2).square()
    cost[n:, m:] = 0
    ii, jj = linear_sum_assignment(cost.detach().cpu().double().numpy())
    sq = cost[ii, jj].sum()
    return sq.clamp_min(1e-20).sqrt() if float(sq.detach()) > 0 else sq


def image_ph(x, y):
    rows = [
        torch.stack([w2(a, b) for a, b in zip(cubical(1 - i), cubical(1 - j))])
        for i, j in zip(x[:, 0], y[:, 0])
    ]
    return torch.stack(rows).mean(0)


def distances(z, normalize=False):
    d = torch.cdist(
        z.flatten(1), z.flatten(1), compute_mode="donot_use_mm_for_euclid_dist"
    )
    if normalize:
        d = d / d.square().mean().clamp_min(1e-12).sqrt()
    return d


def batch_ph(z, received):
    return sum(w2(a, b) for a, b in zip(rips(distances(z)), rips(distances(received))))


def mst_signature(x, z):
    """Bidirectional critical H0 edges, retaining the identities of their endpoints."""
    a, b = distances(x, True), distances(z, True)
    if len(x) < 2:
        return (a.sum() + b.sum()) * 0
    terms = []
    for d in (a, b):
        v = d.detach().cpu().double().numpy()
        tree = minimum_spanning_tree(coo_matrix(np.triu(v + 1.0, k=1)).tocsr()).tocoo()
        terms.append((a[tree.row, tree.col] - b[tree.row, tree.col]).square().mean())
    return sum(terms) / 2
