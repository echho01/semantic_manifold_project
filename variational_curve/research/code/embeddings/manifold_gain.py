"""Held-out point-cloud decoding and nearest-centroid classification.

sigma_total is RMS ambient noise magnitude; per-axis variance is
sigma_total**2 / ambient_dimension. The decoder has reference samples but
never receives the clean held-out source point.
"""
from __future__ import annotations
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import pairwise_distances, pairwise_distances_argmin
from rd import split_indices


def tangent_project(y, X, k=32, pca_dim=2):
    """Project observations onto centred local affine d-planes."""
    y, X = np.asarray(y, float), np.asarray(X, float)
    if not 1 <= pca_dim < min(k, X.shape[1]) or k > len(X):
        raise ValueError("require 1 <= pca_dim < min(k,ambient) and k<=n_reference")
    indices = NearestNeighbors(n_neighbors=k).fit(X).kneighbors(y, return_distance=False)
    out = np.empty_like(y)
    for i, ids in enumerate(indices):
        centre = X[ids].mean(0)
        _, _, vt = np.linalg.svd(X[ids] - centre, full_matrices=False)
        basis = vt[:pca_dim].T
        out[i] = centre + ((y[i] - centre) @ basis) @ basis.T
    return out


def gain_curve(X, levels, n_points=20000, k=32, seed=0, pca_dim=2, projection="tangent"):
    """Operational held-out decoding plus a separately labelled tangent reference."""
    X = np.asarray(X, float)
    train_idx, test_idx = split_indices(len(X), seed)
    reference, source = X[train_idx], X[test_idx[:n_points]]
    n, ambient = source.shape
    k = min(k, len(reference))
    if projection == "tangent":
        decode = lambda y: tangent_project(y, reference, k, pca_dim)
    elif projection == "nearest-sample":
        nearest = NearestNeighbors(n_neighbors=1).fit(reference)
        decode = lambda y: reference[nearest.kneighbors(y, return_distance=False)[:, 0]]
    else:
        raise ValueError("projection must be tangent or nearest-sample")
    zero = float(np.mean(np.sum((decode(source) - source) ** 2, axis=1)))
    out = {"sigma_total": [], "variance_axis": [], "gain_dB": [], "decoded_ambient_mse": [],
           "raw_ambient_mse": [], "zero_noise_decoder_mse": zero,
           "train_count": len(reference), "test_count": len(source), "pca_dim": int(pca_dim),
           "decoder": projection, "oracle_tangent_reference_gain_dB": float(10 * np.log10(ambient / pca_dim)),
           "reference_interpretation": "ideal tangent-noise benchmark; learned affine decoder includes chart and sampling error"}
    rng = np.random.default_rng(seed)
    # Common standardized noise across levels reduces comparison variance.
    noise = rng.normal(size=source.shape)
    for level in levels:
        if level <= 0:
            raise ValueError("noise levels must be positive")
        variance = float(level * level / ambient)
        y = source + noise * np.sqrt(variance)
        decoded = decode(y)
        mse = float(np.mean(np.sum((decoded - source) ** 2, axis=1)))
        raw = float(np.mean(np.sum((y - source) ** 2, axis=1)))
        out["sigma_total"].append(float(level))
        out["variance_axis"].append(variance)
        out["decoded_ambient_mse"].append(mse)
        out["raw_ambient_mse"].append(raw)
        out["gain_dB"].append(float(10 * np.log10(raw / max(mse, 1e-30))))
    return out


def class_separation(X, labels):
    classes = np.unique(labels)
    if len(classes) < 2:
        raise ValueError("at least two classes are required")
    centres = np.stack([X[labels == label].mean(0) for label in classes])
    distances = pairwise_distances(centres)
    return float(distances[np.triu_indices(len(classes), 1)].min())


def outage_curve(X, labels, levels, n_points=20000, seed=0):
    """Train centroids on reference data and classify held-out noisy points."""
    X, labels = np.asarray(X, float), np.asarray(labels)
    if labels.ndim != 1 or len(labels) != len(X):
        raise ValueError("classification labels must align one-to-one with embeddings")
    train_idx, test_idx = split_indices(len(X), seed)
    test_idx = test_idx[:n_points]
    train, train_labels = X[train_idx], labels[train_idx]
    source, truth = X[test_idx], labels[test_idx]
    classes = np.unique(train_labels)
    if len(classes) < 2 or not np.isin(truth, classes).all():
        raise ValueError("all test classes must have training examples")
    centres = np.stack([train[train_labels == label].mean(0) for label in classes])
    def error(y):
        predicted = classes[pairwise_distances_argmin(y, centres)]
        return float(np.mean(predicted != truth))
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=source.shape)
    out = {"sigma_total": [], "variance_axis": [], "classification_error": [],
           "clean_classification_error": error(source),
           "minimum_centroid_separation": class_separation(train, train_labels),
           "train_count": len(train), "test_count": len(source),
           "interpretation": "held-out nearest-centroid classification; centroid separation is not a support-separation bound"}
    for level in levels:
        variance = float(level * level / X.shape[1])
        out["sigma_total"].append(float(level))
        out["variance_axis"].append(variance)
        out["classification_error"].append(error(source + np.sqrt(variance) * noise))
    return out
