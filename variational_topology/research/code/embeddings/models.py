"""Embedding extractors.  All weights download through HuggingFace.

Every extractor exposes `.embed_images(images)` and `.embed_texts(texts)`
and returns float32 arrays of shape [n, D], L2-normalized only when the
model itself does so (CLIP: yes; BERT mean pooling: configurable).
"""
from __future__ import annotations

import os
from typing import Iterable

import numpy as np


def _device():
    import torch
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


class ClipExtractor:
    def __init__(self, hf_id: str, cache_dir: str, batch: int = 64,
                 device: str | None = None):
        import torch
        from transformers import CLIPModel, CLIPProcessor
        self.torch = torch
        self.device = device or _device()
        self.batch = batch
        self.model = CLIPModel.from_pretrained(hf_id, cache_dir=cache_dir).to(self.device).eval()
        self.proc = CLIPProcessor.from_pretrained(hf_id, cache_dir=cache_dir)
        self.D = self.model.config.projection_dim

    def _run(self, fn, inputs):
        import torch
        out = []
        with torch.no_grad():
            for i in range(0, len(inputs), self.batch):
                chunk = inputs[i:i + self.batch]
                feats = fn(**self.proc(text=chunk, return_tensors="pt", padding=True,
                                       truncation=True).to(self.device))
                out.append(feats.float().cpu().numpy())
        return np.concatenate(out, 0)

    def embed_images(self, images: Iterable) -> np.ndarray:
        import torch
        images = list(images)
        out = []
        with torch.no_grad():
            for i in range(0, len(images), self.batch):
                px = self.proc(images=images[i:i + self.batch], return_tensors="pt")["pixel_values"]
                feats = self.model.get_image_features(pixel_values=px.to(self.device))
                out.append(feats.float().cpu().numpy())
        return np.concatenate(out, 0)

    def embed_texts(self, texts: Iterable) -> np.ndarray:
        import torch
        texts = list(texts)
        out = []
        with torch.no_grad():
            for i in range(0, len(texts), self.batch):
                tok = self.proc(text=texts[i:i + self.batch], return_tensors="pt",
                                padding=True, truncation=True, max_length=77)
                feats = self.model.get_text_features(input_ids=tok["input_ids"].to(self.device),
                                                     attention_mask=tok["attention_mask"].to(self.device))
                out.append(feats.float().cpu().numpy())
        return np.concatenate(out, 0)


class BertExtractor:
    def __init__(self, hf_id: str, cache_dir: str, batch: int = 128,
                 device: str | None = None, normalize: bool = False):
        import torch
        from transformers import AutoModel, AutoTokenizer
        self.torch = torch
        self.device = device or _device()
        self.batch = batch
        self.normalize = normalize
        self.tok = AutoTokenizer.from_pretrained(hf_id, cache_dir=cache_dir)
        self.model = AutoModel.from_pretrained(hf_id, cache_dir=cache_dir).to(self.device).eval()
        self.D = self.model.config.hidden_size

    def embed_texts(self, texts: Iterable) -> np.ndarray:
        import torch
        texts = list(texts)
        out = []
        with torch.no_grad():
            for i in range(0, len(texts), self.batch):
                tok = self.tok(texts[i:i + self.batch], return_tensors="pt",
                               padding=True, truncation=True, max_length=256)
                tok = {k: v.to(self.device) for k, v in tok.items()}
                h = self.model(**tok).last_hidden_state            # [B, L, D]
                mask = tok["attention_mask"].unsqueeze(-1).float()
                pooled = (h * mask).sum(1) / mask.sum(1).clamp(min=1e-9)
                if self.normalize:
                    pooled = torch.nn.functional.normalize(pooled, dim=-1)
                out.append(pooled.float().cpu().numpy())
        return np.concatenate(out, 0)

    def embed_images(self, images):
        raise NotImplementedError("BERT is text-only")


class DinoExtractor:
    def __init__(self, hf_id: str, cache_dir: str, batch: int = 64,
                 device: str | None = None):
        import torch
        from transformers import AutoImageProcessor, AutoModel
        self.torch = torch
        self.device = device or _device()
        self.batch = batch
        self.proc = AutoImageProcessor.from_pretrained(hf_id, cache_dir=cache_dir)
        self.model = AutoModel.from_pretrained(hf_id, cache_dir=cache_dir).to(self.device).eval()
        self.D = self.model.config.hidden_size

    def embed_images(self, images) -> np.ndarray:
        import torch
        images = list(images)
        out = []
        with torch.no_grad():
            for i in range(0, len(images), self.batch):
                px = self.proc(images=images[i:i + self.batch], return_tensors="pt")
                h = self.model(**{k: v.to(self.device) for k, v in px.items()})
                out.append(h.last_hidden_state[:, 0].float().cpu().numpy())
        return np.concatenate(out, 0)

    def embed_texts(self, texts):
        raise NotImplementedError("DINOv2 is image-only")


def build(spec: dict, cache_dir: str) -> object:
    kind = spec["kind"]
    if kind == "clip":
        return ClipExtractor(spec["hf_id"], cache_dir)
    if kind == "bert":
        return BertExtractor(spec["hf_id"], cache_dir,
                             normalize=spec.get("normalize", False))
    if kind == "dinov2":
        return DinoExtractor(spec["hf_id"], cache_dir)
    raise ValueError(f"unknown model kind {kind}")
