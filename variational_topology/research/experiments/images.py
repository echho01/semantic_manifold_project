"""Foreground is bright. Splits are made before any resize; no crops cross split boundaries."""

from pathlib import Path
import hashlib
import numpy as np
from PIL import Image, ImageOps
import torch
from torch.utils.data import Dataset


def bucket(name):
    return int.from_bytes(hashlib.sha256(name.encode()).digest()[:4], "big") % 100


class Images(Dataset):
    def __init__(self, root, kind, split, size, limit, seed):
        root = Path(root)
        self.root = root
        self.kind = kind
        self.size = size
        if kind == "omniglot":
            folder = root / (
                "images_evaluation" if split == "test" else "images_background"
            )
            alph = sorted(
                p for p in folder.glob("*") if p.is_dir() and not p.name.startswith(".")
            )
            if split != "test":
                # Reserve whole background alphabets for model selection.
                alph = [
                    p for i, p in enumerate(alph) if (i % 5 == 0) == (split == "val")
                ]
            paths = sorted(p for group in alph for p in group.glob("*/*.png"))
        else:
            # Official challenge test labels are not assumed available. Hash the ORIGINAL labeled tile ID.
            # Supply full-size tiles, named <unique tile ID>_mask.png (or jpg/tif), not pre-cut patches.
            raw = sorted(
                p
                for p in root.rglob("*")
                if p.suffix.lower() in (".png", ".jpg", ".tif", ".tiff")
                and p.stem.lower().endswith("_mask")
            )
            keys = [p.stem[:-5].lower() for p in raw]
            if len(keys) != len(set(keys)):
                raise ValueError(
                    "duplicate road tile IDs; use a single copy of each original labeled tile"
                )
            paths = [
                p
                for key, p in zip(keys, raw)
                if (
                    "train"
                    if bucket(key) < 70
                    else "val"
                    if bucket(key) < 85
                    else "test"
                )
                == split
            ]
        if not paths:
            raise ValueError(f"no {kind} {split} images found below {root}")
        if limit is not None:
            if limit < 1:
                raise ValueError("limit must be positive")
            order = np.random.default_rng(seed).permutation(len(paths))[:limit]
            paths = [paths[i] for i in order]
        self.paths = paths
        self.sha = hashlib.sha256(
            "\n".join(p.relative_to(root).as_posix() for p in paths).encode()
        ).hexdigest()

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, i):
        with Image.open(self.paths[i]) as im:
            im = im.convert("L")
            if self.kind == "omniglot":
                im = ImageOps.invert(im)
            im = im.resize(
                (self.size, self.size),
                Image.Resampling.BILINEAR
                if self.kind == "omniglot"
                else Image.Resampling.NEAREST,
            )
            a = np.array(im, dtype=np.float32) / 255
        return torch.from_numpy(a[None]), self.paths[i].relative_to(
            self.root
        ).as_posix()
