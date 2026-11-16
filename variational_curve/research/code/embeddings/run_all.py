#!/usr/bin/env python3
"""Resumable embedding diagnostics with held-out evaluation and JSON output.

Example: python code/embeddings/run_all.py --config code/embeddings/config.yaml
         --model clip_vit_b32 --dataset cifar10 --stages extract,geometry,rd
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import geometry
import manifold_gain
import rd
import task_mismatch


def load_config(path):
    path = Path(path).resolve()
    if path.suffix in (".yaml", ".yml"):
        import yaml
        cfg = yaml.safe_load(path.read_text())
    else:
        cfg = json.loads(path.read_text())
    for key, value in cfg["paths"].items():
        target = Path(value).expanduser()
        cfg["paths"][key] = str(target if target.is_absolute() else path.parent / target)
    return cfg


def extract(cfg, model_name, dataset_name, out):
    from data_sources import load_dataset_spec, text_for_class_prompt
    from models import build
    model_spec, data_spec = cfg["models"][model_name], cfg["datasets"][dataset_name]
    data = load_dataset_spec(data_spec, cfg["paths"]["torch_cache"], cfg["paths"]["hf_cache"], cfg["seed"])
    model = build(model_spec, cfg["paths"]["hf_cache"])
    emb = {"item_ids": data["item_ids"]}
    modalities = data_spec.get("modalities", ["text"])
    if "image" in modalities and model_spec["kind"] in ("clip", "dinov2"):
        emb["image"] = model.embed_images(data["items"])
    if "text" in modalities and model_spec["kind"] in ("clip", "bert"):
        if data_spec["kind"] in ("cifar10", "cifar100"):
            # One vector per class, not thousands of duplicate pseudo-samples.
            prompts = text_for_class_prompt(data["label_names"], model_spec.get("text_prompt", "a photo of a {label}."))
            emb["text_prototypes"] = model.embed_texts(prompts)
        else:
            emb["text"] = model.embed_texts(data["items"])
    for key in ("labels", "pair_scores", "pair_indices"):
        if data.get(key) is not None:
            emb[key] = data[key]
    if not any(key in emb for key in ("image", "text")):
        raise ValueError("chosen model and dataset have no sample-level common modality")
    np.savez_compressed(out / "arrays.npz", **emb)
    return {"n": len(data["items"]), "label_names": data["label_names"],
            "embedding_shapes": {k: list(v.shape) for k, v in emb.items() if k in ("image", "text", "text_prototypes")},
            "model": model_spec, "dataset": data_spec,
            "resolved_model_revision": getattr(model.model.config, "_commit_hash", None)}


def save_report(out, report):
    temporary = out / "report.json.tmp"
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    temporary.replace(out / "report.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--stages", default="extract,geometry,rd,gain,outage,task")
    parser.add_argument("--model")
    parser.add_argument("--dataset")
    parser.add_argument("--tag")
    parser.add_argument("--modality", choices=("image", "text"))
    args = parser.parse_args()
    cfg = load_config(args.config)
    model_name = args.model or next(iter(cfg["models"]))
    dataset_name = args.dataset or next(iter(cfg["datasets"]))
    tag = args.tag or f"{model_name}__{dataset_name}"
    out = Path(cfg["paths"]["root"]) / tag
    out.mkdir(parents=True, exist_ok=True)
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    unknown = set(stages) - {"extract", "geometry", "rd", "gain", "outage", "task", "mismatch"}
    if unknown:
        raise ValueError(f"unknown stages: {unknown}")
    started = time.time()
    path = out / "report.json"
    report = json.loads(path.read_text()) if path.exists() and "extract" not in stages else {}
    report["config"] = cfg
    report["selection"] = {"model": model_name, "dataset": dataset_name, "tag": tag}
    if "extract" in stages:
        report["extract"] = extract(cfg, model_name, dataset_name, out)
        save_report(out, report)
    array_path = out / "arrays.npz"
    if not array_path.exists():
        raise FileNotFoundError(f"run extract first: {array_path}")
    with np.load(array_path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    key = args.modality or ("image" if "image" in arrays else "text")
    if key not in arrays:
        raise ValueError(f"no sample-level {key} embeddings")
    X = np.asarray(arrays[key], dtype=float)
    labels = arrays.get("labels")
    if not np.isfinite(X).all() or X.ndim != 2:
        raise ValueError("embedding matrix contains nonfinite values or has the wrong shape")
    report["modality"] = key
    report["arrays_sha256"] = hashlib.sha256(array_path.read_bytes()).hexdigest()
    seed = cfg["seed"]
    d_hat = None
    if set(stages) & {"geometry", "rd", "gain"}:
        dim_info = geometry.twonn_dimension(X, seed=seed, **cfg["geometry"].get("twonn", {}))
        d_hat = dim_info["dim"]
        report.setdefault("geometry", {})["twonn"] = dim_info
    if "geometry" in stages:
        params = cfg["geometry"]["curvature"]
        d = params.get("pca_dim") or max(1, int(round(d_hat)))
        if d >= X.shape[1]:
            report["geometry"]["curvature"] = {"status": "estimated dimension fills ambient space; no normal chart"}
        else:
            n_params = 1 + d + d * (d + 1) // 2
            neighbours = max(params.get("k_nn", 120), 2 * n_params)
            neighbours = min(neighbours, len(X) - 1)
            result = geometry.local_curvature(X, d, neighbours, params.get("n_anchor", 4000), seed)
            report["geometry"]["curvature"] = result
            report["geometry"]["continuous_RD_coefficients_nats"] = rd.predict_coefficients(
                result["scal_mean"], d, result["ii2_mean"], result["h2_mean"])
        save_report(out, report)
    if "rd" in stages:
        params = cfg["rd"]
        if params.get("criteria", ["chordal"]) != ["chordal"]:
            raise ValueError("embedding k-means currently evaluates only chordal distortion")
        curve = rd.quantize_curve(X, params["codebook_sizes"], n_fit=params.get("n_fit_points", 200000),
                                  batch=params.get("kmeans_batch", 8192), iters=params.get("kmeans_iters", 50),
                                  seed=seed, test_fraction=params.get("test_fraction", .2))
        report["rd"] = {"curve": curve, "fixed_rate_residual": rd.fit_second_order(curve, d_hat)}
        save_report(out, report)
    if "gain" in stages:
        params = cfg["gain_outage"]
        d = max(1, min(X.shape[1] - 1, int(round(d_hat))))
        report["gain"] = manifold_gain.gain_curve(X, params["noise_levels"], n_points=params.get("n_points", 2000),
                k=max(d + 2, params.get("n_neighbors", 32)), seed=seed, pca_dim=d,
                projection=params.get("projection", "tangent"))
        save_report(out, report)
    if "outage" in stages:
        if labels is None:
            report["outage"] = {"status": "no sample-level classification labels"}
        else:
            params = cfg["gain_outage"]
            report["outage"] = manifold_gain.outage_curve(X, labels, params["noise_levels"],
                                    n_points=params.get("n_points", 2000), seed=seed)
        save_report(out, report)
    if "task" in stages:
        report["task"] = {}
        if labels is not None:
            max_rank = min(X.shape[1], len(np.unique(labels)) - 1)
            for dimension in cfg["task"]["probe_dims"]:
                if dimension > max_rank:
                    report["task"][str(dimension)] = {"status": "requested dimension exceeds class-probe rank"}
                    continue
                curve = task_mismatch.task_rd(X, labels, dimension, cfg["task"]["codebook_sizes"], seed=seed)
                report["task"][str(dimension)] = {"curve": curve, "empirical_exponent": task_mismatch.task_exponent(curve)}
        save_report(out, report)
    if "mismatch" in stages:
        params = cfg["mismatch"]
        other = params["model_b"] if model_name == params["model_a"] else params["model_a"]
        other_path = out.parent / f"{other}__{dataset_name}" / "arrays.npz"
        if not other_path.exists():
            report["mismatch"] = {"status": "other encoder cache absent", "path": str(other_path)}
        else:
            with np.load(other_path, allow_pickle=False) as archive:
                if key not in archive or "item_ids" not in archive or "item_ids" not in arrays:
                    raise ValueError("paired encoder extraction must retain modality and stable item_ids")
                if not np.array_equal(arrays["item_ids"], archive["item_ids"]):
                    raise ValueError("encoder caches have different sample IDs or ordering; re-extract with matching seed/split")
                B = np.asarray(archive[key], float)
            report["mismatch"] = {"paired_edge_deviation": task_mismatch.estimate_epsilon(X, B, seed=seed),
                "curve": task_mismatch.mismatch_rd(X, B, params["codebook_sizes"], seed=seed),
                "alignment": "paired within-encoder distances; orthogonal rotations do not change the metric"}
        save_report(out, report)
    report.setdefault("runs", []).append({"stages": stages, "elapsed_seconds": time.time() - started})
    save_report(out, report)
    print(json.dumps({"report": str(path), "stages": stages, "n": len(X), "ambient_dimension": X.shape[1]}))
    return 0


if __name__ == "__main__":
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())
