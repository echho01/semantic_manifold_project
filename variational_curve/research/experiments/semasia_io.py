"""Pinned SEMASIA Parquet reader; split and model names are part of the sample key."""

import hashlib
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlopen
import numpy as np
import pyarrow.parquet as pq

REPO = "spaicom-lab/semasia-cifar100"
REVISION = "ace0aa91fc0b79953ec1ac619a25072d3a5b3037"
RESNET = "resnet18.a1_in1k"
MOBILE = "mobilenetv3_small_100.lamb_in1k"
HASHES = {
    (
        RESNET,
        "train",
    ): "2dba13a60d7203cba1b1828cb99ce5c64f45c5ebc3624e0ffd2ebd5f9df1c7d9",
    (
        RESNET,
        "test",
    ): "34995557052ba4805960ffdc43b2b218e557a92a086ad77e43c7f38180fb8d36",
    (
        MOBILE,
        "train",
    ): "36f7eb3f307c7be340afe8c5e816ea623c49eb3d37c82e5a067410eb989bf366",
    (
        MOBILE,
        "test",
    ): "0109a7af7d7e3facdd35b0ed7f898f29131f7725bd408bf65ed55ed0f4a8fec2",
}


def stream_hash(stream):
    digest = hashlib.sha256()
    for block in iter(lambda: stream.read(1 << 20), b""):
        digest.update(block)
    return digest.hexdigest()


def shard_file(directory, model, split, download):
    relative = Path(split) / model / "part-00000.parquet"
    path = directory / relative
    if not path.exists() and download:
        path.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://huggingface.co/datasets/{REPO}/resolve/{REVISION}/{relative.as_posix()}?download=true"
        with (
            urlopen(url, timeout=90) as response,
            NamedTemporaryFile(dir=path.parent, delete=False) as tmp,
        ):
            temp_path = Path(tmp.name)
            try:
                while block := response.read(1024 * 1024):
                    tmp.write(block)
                tmp.close()
                with temp_path.open("rb") as stream:
                    if stream_hash(stream) != HASHES[model, split]:
                        raise ValueError(f"unexpected download for {relative}")
                temp_path.replace(path)
            finally:
                tmp.close()
                temp_path.unlink(missing_ok=True)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Use --download to fetch the pinned files."
        )
    with path.open("rb") as stream:
        if stream_hash(stream) != HASHES[model, split]:
            raise ValueError(f"file differs from the pinned revision: {path}")
    return pq.ParquetFile(path)


def metadata(shard, model):
    t = shard.read(columns=["id", "fine_label", "coarse_label", "model_name"])
    if set(t["model_name"].unique().to_pylist()) != {model}:
        raise ValueError("unexpected model_name in shard")
    ids, fine, coarse = (
        t[name].to_numpy() for name in ("id", "fine_label", "coarse_label")
    )
    order = np.argsort(ids)
    if len(np.unique(ids)) != len(ids):
        raise ValueError("duplicate sample IDs in one split")
    return np.column_stack((ids[order], fine[order], coarse[order]))


def features(shard, ids):
    """Reorder selected rows by sample ID, even if different models store different row orders."""
    ids = np.asarray(ids)
    sorted_rows = np.argsort(ids)
    sorted_ids = ids[sorted_rows]
    d = shard.schema_arrow.field("embedding").type.list_size
    out = np.empty((len(ids), d), dtype=np.float64)
    seen = np.zeros(len(ids), dtype=bool)
    for batch in shard.iter_batches(batch_size=2048, columns=["id", "embedding"]):
        batch_ids = batch.column(0).to_numpy()
        slots = np.minimum(np.searchsorted(sorted_ids, batch_ids), len(ids) - 1)
        keep = sorted_ids[slots] == batch_ids
        rows = sorted_rows[slots[keep]]
        vectors = batch.column(1).flatten().to_numpy().reshape(len(batch_ids), d)
        out[rows] = vectors[keep]
        seen[rows] = True
    if not seen.all() or not np.isfinite(out).all():
        raise ValueError("missing or nonfinite feature")
    norms = np.linalg.norm(out, axis=1, keepdims=True)
    if np.any(norms == 0):
        raise ValueError("cosine geometry requires nonzero vectors")
    return out / norms
