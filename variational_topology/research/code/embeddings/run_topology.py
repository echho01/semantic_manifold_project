#!/usr/bin/env python
"""Topology measurements on paired cached embeddings with held-out receivers.

Caches must contain item_ids. Rows are aligned by identity, never by their
position or the length of the shortest cache. Class means, normalization,
and the scalar probe use a disjoint reference split. No model is loaded.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.model_selection import train_test_split
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).parent))
import topology_exp as T


def load_config(path: str) -> dict:
    path = Path(path)
    if path.suffix in {".yaml", ".yml"}:
        import yaml
        return yaml.safe_load(path.read_text())
    return json.loads(path.read_text())


def load(root: Path, model: str, dataset: str, key: str | None):
    path = root / f"{model}__{dataset}" / "arrays.npz"
    with np.load(path, allow_pickle=False) as z:
        k = key or ("image" if "image" in z.files else "text")
        if "item_ids" not in z.files:
            raise ValueError(f"{path}: missing item_ids; rerun extraction to establish row identity")
        X = np.asarray(z[k], dtype=np.float64)
        ids = np.asarray(z["item_ids"]).astype(str)
        labels = np.asarray(z["labels"]) if "labels" in z.files else None
    if X.ndim != 2 or len(ids) != len(X) or len(np.unique(ids)) != len(ids):
        raise ValueError(f"{path}: invalid feature rows or duplicate item_ids")
    if not np.isfinite(X).all() or (labels is not None and len(labels) != len(X)):
        raise ValueError(f"{path}: invalid feature values or labels")
    return {"X": X, "ids": ids, "labels": labels,
            "path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def align(raw: dict):
    names = list(raw)
    common = set(raw[names[0]]["ids"])
    for r in raw.values():
        common.intersection_update(r["ids"])
    ids = np.asarray(sorted(common))
    if len(ids) < 4:
        raise ValueError("fewer than four shared item identities")
    aligned, labels = {}, None
    for name, r in raw.items():
        lookup = {s: i for i, s in enumerate(r["ids"])}
        ix = np.asarray([lookup[s] for s in ids])
        aligned[name] = r["X"][ix]
        y = r["labels"][ix] if r["labels"] is not None else None
        if y is not None:
            if labels is not None and not np.array_equal(y, labels):
                raise ValueError(f"labels disagree across caches after ID alignment: {name}")
            labels = y
    return ids, aligned, labels


def run(cfg: dict, dataset: str, models: list, stages: list, modality=None):
    tc = cfg["topology"]
    root = Path(cfg["paths"]["root"])
    seed = int(cfg["seed"])
    t0 = time.time()
    raw = {m: load(root, m, dataset, modality) for m in models}
    ids, aligned, labels = align(raw)
    stratify = labels if labels is not None else None
    reference, evaluation = train_test_split(np.arange(len(ids)),
        test_size=tc.get("test_fraction", 0.4), random_state=seed, stratify=stratify)
    reference, evaluation = np.sort(reference), np.sort(evaluation)
    items = evaluation[T.common_items(len(evaluation), tc["n_items"], seed)]
    train, heldout, clouds, scales = {}, {}, {}, {}
    for m, X in aligned.items():
        mean, scale = T.normalization(X[reference], seed)
        train[m] = (X[reference] - mean) / scale
        heldout[m] = (X[evaluation] - mean) / scale
        clouds[m] = (X[items] - mean) / scale
        scales[m] = float(scale)
    yref = labels[reference] if labels is not None else None
    yeval = labels[evaluation] if labels is not None else None
    report = {"config": {"dataset": dataset, "models": models, "stages": stages,
        "seed": seed, "n_items": len(items), "D_a": {m: clouds[m].shape[1] for m in models},
        "topology": tc}, "provenance": {
        "inputs": {m: {k: r[k] for k in ["path", "sha256"]} for m, r in raw.items()},
        "reference_ids": ids[reference].tolist(), "evaluation_ids": ids[evaluation].tolist(),
        "diagram_ids": ids[items].tolist(), "normalization_reference_median": scales,
        "item_alignment": "intersection_by_unique_ID_then_label_agreement",
        "source_sha256": {str(p.name): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__), Path(T.__file__)]}}}
    if "floor" in stages and len(models) >= 2:
        report["floor"] = T.floor_matrix(clouds, tc["maxdim"])
    if "dictionary" in stages and len(models) >= 2:
        sizes = [s for s in tc["dictionary_sizes"] if s <= len(items)]
        if not sizes:
            raise ValueError("all dictionary sizes exceed the held-out cloud size")
        report["dictionary"] = {f"{a}->{b}": T.dictionary_curve(
            clouds[a], clouds[b], sizes, tc["maxdim"], seed)
            for a, b in itertools.permutations(models, 2)}
    if "reparam" in stages:
        report["reparam"] = {m: T.reparam_check(clouds[m][:tc["reparam_points"]],
            tc["noise_levels"], tc["maxdim"], seed=seed) for m in models}
    if "cliff" in stages:
        if labels is None:
            raise ValueError("cliff stage requires aligned class labels")
        report["cliff"] = {m: T.cliff_curve(heldout[m], yeval,
            tc["cliff_noise_levels"], per_class=tc["cliff_per_class"], seed=seed,
            X_reference=train[m], labels_reference=yref) for m in models}
    if "reeb" in stages:
        if labels is None:
            raise ValueError("reeb probe stage requires aligned class labels")
        report["reeb"] = {}
        for m in models:
            f = T.probe_score(clouds[m], X_reference=train[m], labels_reference=yref)
            report["reeb"][m] = T.reeb_overhead(clouds[m], f, n_bins=tc["reeb_bins"])
    report["elapsed_s"] = time.time() - t0
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--models", required=True)
    ap.add_argument("--modality")
    ap.add_argument("--stages", default="floor,dictionary,reparam,cliff,reeb")
    ap.add_argument("--tag")
    args = ap.parse_args()
    cfg = load_config(args.config)
    stages = args.stages.split(",")
    if set(stages) - {"floor", "dictionary", "reparam", "cliff", "reeb"}:
        ap.error("unknown topology stage")
    with threadpool_limits(limits=1):
        report = run(cfg, args.dataset, args.models.split(","), stages, args.modality)
    out = Path(cfg["paths"]["root"]) / (args.tag or f"topology__{args.dataset}")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Aligned {report['config']['n_items']} held-out items; report: {out / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
