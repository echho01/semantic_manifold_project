"""Paired CIFAR-100 feature pilot with an actual fixed-length index stream.

The input directory is the official split/model/part-00000.parquet layout.
Use --download to fetch only the four pinned files into that directory.
"""
import argparse
import hashlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlopen

import numpy as np
import pyarrow.parquet as pq
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp, softmax
from sklearn.linear_model import LogisticRegression
from threadpoolctl import threadpool_limits

from manifold_rd import finite_ba, posterior_kl_bits

ROOT = Path(__file__).resolve().parent
from geometry import twonn_dimension

REPO = "spaicom-lab/semasia-cifar100"
REVISION = "ace0aa91fc0b79953ec1ac619a25072d3a5b3037"
RESNET = "resnet18.a1_in1k"
MOBILE = "mobilenetv3_small_100.lamb_in1k"
HASHES = {
    (RESNET, "train"): "2dba13a60d7203cba1b1828cb99ce5c64f45c5ebc3624e0ffd2ebd5f9df1c7d9",
    (RESNET, "test"): "34995557052ba4805960ffdc43b2b218e557a92a086ad77e43c7f38180fb8d36",
    (MOBILE, "train"): "36f7eb3f307c7be340afe8c5e816ea623c49eb3d37c82e5a067410eb989bf366",
    (MOBILE, "test"): "0109a7af7d7e3facdd35b0ed7f898f29131f7725bd408bf65ed55ed0f4a8fec2",
}
SIZES = (1, 4, 16, 64, 256)


def sha256_file(path):
    """Streamed digest; hashlib.file_digest needs Python 3.11."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def shard_file(directory, model, split, download):
    relative = Path(split) / model / "part-00000.parquet"
    path = directory / relative
    if not path.exists() and download:
        path.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://huggingface.co/datasets/{REPO}/resolve/{REVISION}/{relative.as_posix()}?download=true"
        with urlopen(url, timeout=90) as response, NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
            temp_path = Path(tmp.name)
            try:
                while block := response.read(1024 * 1024):
                    tmp.write(block)
                tmp.close()
                if sha256_file(temp_path) != HASHES[model, split]:
                    raise ValueError(f"unexpected download for {relative}")
                temp_path.replace(path)
            finally:
                tmp.close()
                temp_path.unlink(missing_ok=True)
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Use --download to fetch the pinned files.")
    if sha256_file(path) != HASHES[model, split]:
        raise ValueError(f"file differs from the pinned revision: {path}")
    return pq.ParquetFile(path)


def metadata(shard, model):
    t = shard.read(columns=["id", "fine_label", "coarse_label", "model_name"])
    if set(t["model_name"].unique().to_pylist()) != {model}:
        raise ValueError("unexpected model_name in shard")
    ids, fine, coarse = (t[name].to_numpy() for name in ("id", "fine_label", "coarse_label"))
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


def task_readout(train, labels, n_fit, test):
    head = LogisticRegression(C=10.0, max_iter=600, tol=1e-7).fit(train[:n_fit], labels[:n_fit])
    if not np.array_equal(head.classes_, np.arange(20)):
        raise ValueError("the fitting subset must include all 20 coarse classes")
    logits = head.decision_function(train[n_fit:])
    target = labels[n_fit:]

    def loss(log_temp):
        scores = logits / np.exp(log_temp)
        return float(np.mean(logsumexp(scores, axis=1) - scores[np.arange(len(target)), target]))

    fitted = minimize_scalar(loss, bounds=(np.log(0.1), np.log(10)), method="bounded")
    if not fitted.success:
        raise RuntimeError("temperature fit failed")
    temp = float(np.exp(fitted.x))
    q_train = softmax(head.decision_function(train[:n_fit]) / temp, axis=1)
    q_test = softmax(head.decision_function(test) / temp, axis=1)
    return q_train, q_test, temp


def label_scores(q, labels):
    return dict(accuracy=float(np.mean(q.argmax(1) == labels)),
                log_loss_bits=float(-np.log2(q[np.arange(len(labels)), labels]).mean()))


def transmit_indices(indices, size):
    """Fixed-width binary indices, with the last byte padded if necessary."""
    if size < 1 or np.any(indices < 0) or np.any(indices >= size):
        raise ValueError("index outside the message alphabet")
    width = (size - 1).bit_length()
    if width == 0:
        return np.zeros_like(indices), 0
    shifts = np.arange(width - 1, -1, -1)
    bits = ((indices[:, None] >> shifts) & 1).astype(np.uint8)
    payload = np.packbits(bits)
    decoded = np.unpackbits(payload, count=len(indices) * width).reshape(-1, width) @ (1 << shifts)
    np.testing.assert_array_equal(decoded, indices)
    return decoded, int(payload.nbytes)


def run(directory, download, seed):
    rng = np.random.default_rng(seed)
    n_fit, n_cal, n_test = 8192, 2048, 2048
    report = dict(dataset=REPO, revision=REVISION, seed=seed, fit=n_fit, calibration=n_cal,
                  test=n_test, task="coarse_label", models={}, links=[], task_baselines=[])
    reference, selection, clouds = {}, {}, {}
    for model in (RESNET, MOBILE):
        arrays, labels, row_counts = {}, {}, {}
        for split, count in (("train", n_fit + n_cal), ("test", n_test)):
            shard = shard_file(directory, model, split, download)
            meta = metadata(shard, model)
            if split not in reference:
                reference[split] = meta
                selection[split] = rng.choice(len(meta), count, replace=False)
            np.testing.assert_array_equal(meta, reference[split])
            selected = meta[selection[split]]
            arrays[split] = features(shard, selected[:, 0])
            labels[split] = selected[:, 2]
            row_counts[split] = len(meta)
        train, test = arrays["train"], arrays["test"]
        q_train, q_test, temp = task_readout(train, labels["train"], n_fit, test)
        book = train[:max(SIZES)]
        book_q = q_train[:max(SIZES)]
        costs = posterior_kl_bits(q_train[256:384], book_q[:16])
        ba_points = []
        for beta in (0, 0.5, 1, 2, 4):
            a = finite_ba(np.full(len(costs), 1 / len(costs)), costs, beta)
            ba_points.append(dict(beta=beta, rate_bits=a.rate, posterior_kl_bits=a.distortion,
                                  lagrangian_gap_bits=a.gap))
        clouds[model] = (test, q_test, book, book_q)
        report["models"][model] = dict(
            dimension=test.shape[1], shard_rows=row_counts, temperature=temp,
            clean=label_scores(q_test, labels["test"]), finite_ba_128_by_16=ba_points,
            twonn={n: twonn_dimension(test[:n]) for n in (512, 2048)})
    report["id_sha256"] = {}
    for split, selected in selection.items():
        ids = reference[split][selected, 0].astype("<u4")
        report["id_sha256"][split] = hashlib.sha256(ids.tobytes()).hexdigest()
    test_labels = reference["test"][selection["test"], 2]
    for sender, (src, sender_q, anchors, sender_book_q) in clouds.items():
        similarity = src @ anchors.T
        task_cost = posterior_kl_bits(sender_q, sender_book_q)
        for size in SIZES:
            received, n_bytes = transmit_indices(similarity[:, :size].argmax(1), size)
            for receiver, (dst, q, dst_book, dst_book_q) in clouds.items():
                pred = dst_book_q[received]
                kl = np.sum(q * (np.log2(q) - np.log2(pred)), axis=1)
                report["links"].append(dict(
                    sender=sender, receiver=receiver, alphabet=size, payload_bytes=n_bytes,
                    bits_per_image=8*n_bytes/len(src), mean_posterior_kl_bits=float(kl.mean()),
                    mean_squared_chord=float(np.mean(np.sum((dst - dst_book[received])**2, axis=1))),
                    **label_scores(pred, test_labels)))
            task_index, task_bytes = transmit_indices(task_cost[:, :size].argmin(1), size)
            report["task_baselines"].append(dict(
                sender=sender, method="same_model_KL_nearest_anchor", alphabet=size,
                bits_per_image=8*task_bytes/len(src),
                mean_posterior_kl_bits=float(np.mean(task_cost[np.arange(len(src)), task_index])),
                **label_scores(sender_book_q[task_index], test_labels)))
        decisions, decision_bytes = transmit_indices(sender_q.argmax(1), sender_q.shape[1])
        report["task_baselines"].append(dict(
            sender=sender, method="direct_coarse_decision", alphabet=sender_q.shape[1],
            bits_per_image=8*decision_bytes/len(src), accuracy=float(np.mean(decisions == test_labels))))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="local cache of the pinned Parquet files")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "semasia_cifar100.json")
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        result = run(args.data, args.download, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Matched both official splits. Wrote {len(result['links'])} codec measurements to {args.output}")
