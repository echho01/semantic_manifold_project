#!/usr/bin/env python
"""Paper 1 experiment driver.

Usage:
    python run_all.py --config config.yaml
    python run_all.py --config config.yaml --stages geometry,rd
    python run_all.py --config config.yaml --model clip_vit_b32 --dataset cifar10

Stages: extract, geometry, rd, gain, outage, task, mismatch.
Outputs: <root>/<tag>/<stage>.json and <root>/<tag>/arrays.npz
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import geometry
import manifold_gain
import rd
import task_mismatch
from data_sources import load_dataset_spec
from models import build

try:
    import yaml
except ImportError:      # allow running without PyYAML by using JSON configs
    yaml = None


def load_config(path: str) -> dict:
    text = Path(path).read_text()
    if path.endswith((".yaml", ".yml")):
        if yaml is None:
            raise SystemExit("PyYAML not installed; use a .json config or pip install pyyaml")
        return yaml.safe_load(text)
    return json.loads(text)


def extract(cfg: dict, model_name: str, dataset_name: str, out: Path) -> dict:
    ma = cfg["models"][model_name]
    da = cfg["datasets"][dataset_name]
    cache = cfg["paths"]["hf_cache"]
    ds = load_dataset_spec(da, cfg["paths"]["torch_cache"], cache, seed=cfg["seed"])
    model = build(ma, cache)
    emb = {}
    if "image" in da.get("modalities", []) and hasattr(model, "embed_images"):
        try:
            emb["image"] = model.embed_images(ds["items"])
        except NotImplementedError:
            pass
    if "text" in da.get("modalities", []) and hasattr(model, "embed_texts"):
        try:
            if ma.get("text_prompt") and ds.get("label_names"):
                from data_sources import text_for_class_prompt
                texts = [text_for_class_prompt(ds["label_names"], ma["text_prompt"])[
                    int(y)] for y in ds["labels"]]
            else:
                texts = ds["items"]
            emb["text"] = model.embed_texts(texts)
        except NotImplementedError:
            pass
    emb["item_ids"] = np.asarray(ds["item_ids"])
    emb["labels"] = ds["labels"]
    emb["label_names"] = ds["label_names"]
    np.savez_compressed(out / "arrays.npz", **{k: v for k, v in emb.items()
                                               if isinstance(v, np.ndarray)})
    return {"n": len(ds["items"]), "modalities": [k for k in emb if isinstance(emb[k], np.ndarray)],
            "D": {k: int(emb[k].shape[1]) for k in emb if isinstance(emb[k], np.ndarray)
                  and getattr(emb[k], "ndim", 0) == 2}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--stages", default="extract,geometry,rd,gain,outage,task,mismatch")
    ap.add_argument("--model", default=None)
    ap.add_argument("--dataset", default=None)
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    model_name = args.model or list(cfg["models"])[0]
    dataset_name = args.dataset or list(cfg["datasets"])[0]
    tag = args.tag or f"{model_name}__{dataset_name}"
    out = Path(cfg["paths"]["root"]) / tag
    out.mkdir(parents=True, exist_ok=True)
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]

    t0 = time.time()
    report = {"config": {"model": model_name, "dataset": dataset_name,
                         "tag": tag, "stages": stages}}

    if "extract" in stages:
        report["extract"] = extract(cfg, model_name, dataset_name, out)

    Z = np.load(out / "arrays.npz", allow_pickle=True) if (out / "arrays.npz").exists() else None
    if Z is None:
        print("no embeddings found; run the extract stage first", file=sys.stderr)
        return 1

    key = "image" if "image" in Z.files else ("text" if "text" in Z.files else None)
    if key is None:
        print("arrays.npz has no embedding matrix", file=sys.stderr)
        return 1
    X = Z[key].astype(np.float32)
    labels = Z["labels"] if "labels" in Z.files else None

    # geometry
    gcfg = cfg["geometry"]
    dim_info = geometry.twonn_dimension(X, **{k: v for k, v in gcfg["twonn"].items()
                                              if k in ("max_points", "frac")},
                                        seed=cfg["seed"])
    d_hat = dim_info["dim"]
    report["geometry"] = {"twonn": dim_info}

    if "geometry" in stages or "rd" in stages or "gain" in stages:
        ccfg = dict(gcfg["curvature"])
        if ccfg.get("pca_dim") is None:
            ccfg["pca_dim"] = max(1, int(round(d_hat)))
        curv = geometry.local_curvature(X, **{k: v for k, v in ccfg.items()
                                              if k in ("pca_dim", "k_nn", "n_anchor")},
                                        seed=cfg["seed"])
        report["geometry"]["curvature"] = curv
        tau_hat = curv["reach_median"]

    # rate-distortion
    if "rd" in stages:
        rcfg = cfg["rd"]
        curves = {}
        for crit in rcfg["criteria"]:
            c = rd.quantize_curve(X, rcfg["codebook_sizes"], n_fit=rcfg["n_fit_points"],
                                  batch=rcfg["kmeans_batch"], iters=rcfg["kmeans_iters"],
                                  seed=cfg["seed"], geodesic=(crit == "geodesic"),
                                  tau=tau_hat if rcfg.get("arc_correction") else None)
            fit = rd.fit_second_order(c, d_hat)
            curves[crit] = {"curve": {k: v for k, v in c.items() if k != "labels"},
                            "fit": fit}
        pred = rd.predict_coefficients(curv["scal_mean"], d_hat,
                                       anisotropy=curv["ii2_mean"])
        report["rd"] = {"curves": curves, "prediction_nats": pred}

    # gain
    if "gain" in stages:
        gcfg2 = cfg["gain_outage"]
        report["gain"] = manifold_gain.gain_curve(
            X, gcfg2["noise_levels"], n_points=gcfg2["n_points"],
            k=gcfg2["n_neighbors"], seed=cfg["seed"])
        report["gain"]["theory_dB"] = 10 * np.log10(X.shape[1] / d_hat)
        report["gain"]["d_hat"] = d_hat

    # outage
    if "outage" in stages and labels is not None:
        gcfg2 = cfg["gain_outage"]
        report["outage"] = manifold_gain.outage_curve(
            X, labels, gcfg2["noise_levels"], n_points=gcfg2["n_points"],
            seed=cfg["seed"])
        report["outage"]["reach_hat"] = tau_hat

    # task
    if "task" in stages and labels is not None:
        tcfg = cfg["task"]
        report["task"] = {}
        for dprobe in tcfg["probe_dims"]:
            if dprobe > X.shape[1]:
                continue
            c = task_mismatch.task_rd(X, labels, dprobe, tcfg["codebook_sizes"],
                                      seed=cfg["seed"])
            report["task"][str(dprobe)] = {"curve": c,
                                           "exponent": task_mismatch.task_exponent(c)}

    # mismatch
    if "mismatch" in stages:
        mcfg = cfg["mismatch"]
        other = mcfg["model_b"] if mcfg["model_a"] == model_name else mcfg["model_a"]
        rep_path = out.parent / f"{other}__{dataset_name}" / "arrays.npz"
        if rep_path.exists():
            Z2 = np.load(rep_path, allow_pickle=True)
            if key in Z2.files:
                B = Z2[key].astype(np.float32)
                n = min(len(X), len(B))
                A_, B_ = X[:n], B[:n]
                if mcfg["align"] == "procrustes":
                    # Procrustes needs equal dimensions: project both to the
                    # smaller one by PCA
                    k = min(A_.shape[1], B_.shape[1])
                    A_ = _pca(A_, k); B_ = _pca(B_, k)
                    B_ = task_mismatch.procrustes(A_, B_)
                eps = task_mismatch.estimate_epsilon(A_, B_, seed=cfg["seed"])
                mm = task_mismatch.mismatch_rd(A_, B_, mcfg["codebook_sizes"],
                                               seed=cfg["seed"])
                report["mismatch"] = {"epsilon_g": eps, "curve": mm,
                                      "prediction": "ratio - 1 <= C * eps_g, flat in D"}
        else:
            report["mismatch"] = {"skipped": f"{rep_path} not found"}

    report["elapsed_s"] = time.time() - t0
    (out / "report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "rd"}, indent=2)[:4000])
    print(f"\nfull report: {out / 'report.json'}")
    return 0


def _pca(X: np.ndarray, k: int) -> np.ndarray:
    mu = X.mean(0)
    Xc = X - mu
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return (U[:, :k] * S[:k])


if __name__ == "__main__":
    raise SystemExit(main())
