"""Dataset adapters with stable item IDs and separate pair-score metadata.

This module is deliberately not named datasets.py: that name shadows the
Hugging Face datasets package when scripts run from this directory.
"""
from __future__ import annotations
import numpy as np


def _cifar(kind, split, max_items, torch_cache, seed=0):
    from torchvision.datasets import CIFAR10, CIFAR100
    cls = CIFAR10 if kind == "cifar10" else CIFAR100
    data = cls(root=torch_cache, train=(split == "train"), download=True, transform=None)
    ids = np.random.default_rng(seed).choice(len(data), min(max_items, len(data)), replace=False)
    items, labels = [], []
    for index in ids:
        image, label = data[int(index)]
        items.append(image.convert("RGB"))
        labels.append(int(label))
    return {"items": items, "labels": np.asarray(labels, dtype=np.int64),
            "label_names": list(data.classes), "item_ids": ids.astype(np.int64)}


def _hf_text(spec, hf_cache, seed=0):
    from datasets import load_dataset
    data = load_dataset(spec["hf_id"], split=spec["split"], cache_dir=hf_cache,
                        revision=spec.get("revision"))
    ids = np.random.default_rng(seed).choice(len(data), min(spec.get("max_items", len(data)), len(data)), replace=False)
    data = data.select(ids.tolist())
    field = spec.get("label_field")
    labels = np.asarray(data[field]) if field in data.column_names else None
    feature = data.features[field] if labels is not None else None
    names = getattr(feature, "names", None)
    return {"items": [str(t) for t in data[spec["text_field"]]], "labels": labels,
            "label_names": names, "item_ids": ids.astype(np.int64)}


def _hf_pairs(spec, hf_cache, seed=0):
    from datasets import load_dataset
    data = load_dataset(spec["hf_id"], split=spec["split"], cache_dir=hf_cache,
                        revision=spec.get("revision"))
    ids = np.random.default_rng(seed).choice(len(data), min(spec.get("max_items", len(data)), len(data)), replace=False)
    data = data.select(ids.tolist())
    a = [str(t) for t in data[spec.get("first_field", "sentence1")]]
    b = [str(t) for t in data[spec.get("second_field", "sentence2")]]
    return {"items": a + b, "labels": None, "label_names": None,
            "item_ids": np.asarray([f"a:{i}" for i in ids] + [f"b:{i}" for i in ids]),
            "pair_scores": np.asarray(data[spec.get("score_field", "score")], dtype=float),
            "pair_indices": np.column_stack([np.arange(len(a)), np.arange(len(a)) + len(a)])}


def load_dataset_spec(spec, torch_cache, hf_cache, seed=0):
    kind = spec["kind"]
    if kind in ("cifar10", "cifar100"):
        return _cifar(kind, spec.get("split", "test"), spec["max_items"], torch_cache, seed)
    if kind == "hf_text":
        return _hf_text(spec, hf_cache, seed)
    if kind == "hf_pairs":
        return _hf_pairs(spec, hf_cache, seed)
    raise ValueError(f"unknown dataset kind {kind}")


def text_for_class_prompt(label_names, prompt):
    return [prompt.format(label=name.replace("_", " ")) for name in label_names]
