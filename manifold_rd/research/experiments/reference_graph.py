"""Training-only graph approximation to a fixed task geometry, with out-of-sample attachment."""

import math
import numpy as np
import torch
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra


def coordinates(x, q, eta):
    return torch.cat((2 * q.sqrt(), math.sqrt(eta) * x), 1)


def neighbors(x, reference, k, exclude_self=False):
    vals = []
    idx = []
    for start in range(0, len(x), 256):
        block = x[start : start + 256]
        dist = torch.cdist(
            block, reference, compute_mode="donot_use_mm_for_euclid_dist"
        )
        if exclude_self:
            v = torch.arange(len(block), device=x.device)
            dist[v, start + v] = float("inf")
        val, ix = dist.topk(k, largest=False)
        vals.append(val)
        idx.append(ix)
    return torch.cat(vals), torch.cat(idx)


def fit_graph(z, book, k):
    if not 1 <= k < len(z):
        raise ValueError("graph k must be in [1,fit-1]")
    d, ix = neighbors(z, z, k, True)
    ii = np.repeat(np.arange(len(z)), k)
    jj = ix.cpu().numpy().ravel()
    # Weighted union kNN graph: build exactly one copy of each undirected edge.
    edge = np.column_stack((np.minimum(ii, jj), np.maximum(ii, jj)))
    edge, loc = np.unique(edge, axis=0, return_index=True)
    w = d.cpu().double().numpy().ravel()[loc]
    a, b = edge.T
    graph = csr_matrix(
        (np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(len(z), len(z))
    )
    n, labels = connected_components(graph, directed=False)
    if n != 1:
        raise ValueError(
            f"training graph has {n} components; choose a documented larger k or an explicit component model"
        )
    length = np.ascontiguousarray(
        dijkstra(graph, directed=False, indices=np.arange(book)).T
    )
    table = torch.from_numpy(length)
    return dict(
        points=z.cpu(),
        to_anchor=table,
        k=k,
        edges=len(w),
        convention="union kNN, Euclidean length in frozen task coordinates, squared shortest path",
    ), table.to(z.device).square()


def attached(z, graph):
    pts = graph["points"].to(z.device)
    table = graph["to_anchor"].to(z.device)
    out = []
    for batch in z.split(128):
        d, i = neighbors(batch, pts, graph["k"])
        out.append((d.double()[:, :, None] + table[i]).amin(1).square())
    return torch.cat(out)
