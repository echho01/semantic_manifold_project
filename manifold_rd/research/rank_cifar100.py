"""Task-rank law on real posterior clouds.

By posterior sufficiency the task rate-distortion function of a representation U is that of
its posterior Q = q(U) on the simplex. A c-way task has rank r <= c - 1, so the law predicts
that the rate grows by r/2 bits for each halving of the KL distortion. A cloud of a few
thousand points can resolve that slope when r is 1, 2 or 3, which it cannot for the
representation itself.

Tasks are built by merging the 20 CIFAR-100 superclasses into c groups. For each task a
calibrated c-way head on a frozen backbone gives the posteriors of official test images. The
source is the equiprobable cloud of n posteriors and the object receiver returns one of them.
Every point is an achievable Blahut point with its dual value. An empirical cloud of n points
cannot carry more than log2(n) bits, so the slope is read only where the clouds of n and n/2
points agree.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from knowledge_cifar100 import calibrated_head
from manifold_rd import finite_ba, posterior_kl_bits
from semasia_cifar100 import MOBILE, REPO, RESNET, REVISION, features, metadata, shard_file

ROOT = Path(__file__).resolve().parent


def ba_torch(cost, beta, device, tol=1e-5, max_iter=3000):
    """The iteration of manifold_rd.finite_ba for a uniform source, on a GPU."""
    import torch
    c = torch.as_tensor(cost, dtype=torch.float64, device=device)
    row_min = c.min(1, keepdim=True).values
    kernel = torch.exp2(-beta * (c - row_min))
    prior = torch.full((c.shape[1],), 1.0 / c.shape[1], dtype=torch.float64, device=device)
    for _ in range(max_iter):
        z = kernel @ prior
        gain = (kernel / z[:, None]).mean(0)
        log_gain = torch.log2(gain.max())
        if log_gain <= tol:
            break
        prior = (prior * gain).clamp_min(1e-250)
        prior = prior / prior.sum()
    w = kernel * prior / z[:, None]
    marginal = w.mean(0)
    rate = (w * (torch.log2(w.clamp_min(1e-300)) - torch.log2(marginal.clamp_min(1e-300)))).sum(1).mean()
    distortion = (w * c).sum(1).mean()
    dual = beta * row_min.mean() - torch.log2(z).mean() - log_gain
    return float(rate), float(distortion), max(0.0, float(rate + beta * distortion - dual))


def curve(q, betas, device=None):
    cost = posterior_kl_bits(q, q)
    if device:
        return [dict(zip(("beta", "rate_bits", "kl_bits", "gap_bits"), (float(b), *ba_torch(cost, float(b), device))))
                for b in betas]
    mass = np.full(len(q), 1.0 / len(q))
    points = [finite_ba(mass, cost, b, tol=1e-5, max_iter=3000) for b in betas]
    return [dict(beta=float(a.beta), rate_bits=a.rate, kl_bits=a.distortion, gap_bits=a.gap) for a in points]


def rate_at(points, distortion):
    d = np.array([p["kl_bits"] for p in points])[::-1]
    r = np.array([p["rate_bits"] for p in points])[::-1]
    return float(np.interp(np.log(distortion), np.log(d), r))


def slope_summary(full, half, tolerance):
    """Bits per halving of the distortion over the range where n and n/2 points agree."""
    lo = max(full[-1]["kl_bits"], half[-1]["kl_bits"])
    hi = min(full[0]["kl_bits"], half[0]["kl_bits"])
    grid = np.geomspace(hi, lo, 200)
    r_full = np.array([rate_at(full, d) for d in grid])
    r_half = np.array([rate_at(half, d) for d in grid])
    ok = (np.abs(r_full - r_half) <= tolerance) & (r_full >= 1.0)
    if ok.sum() < 5:
        return None
    # Longest run of agreement starting from the low-rate end.
    last = np.flatnonzero(ok)[-1] if ok[np.flatnonzero(ok)[0]:].all() else np.flatnonzero(~ok & (r_full >= 1.0))[0] - 1
    first = np.flatnonzero(ok)[0]
    if last - first < 4:
        return None
    d1, d2 = grid[first], grid[last]
    r1, r2 = r_full[first], r_full[last]
    mid = (first + last) // 2  # The finer half of the range is closer to the high-resolution regime.
    return dict(distortion_range=[float(d1), float(d2)], rate_range=[float(r1), float(r2)],
                bits_per_halving=float((r2 - r1) / np.log2(d1 / d2)),
                bits_per_halving_finer_half=float((r2 - r_full[mid]) / np.log2(grid[mid] / d2)))


def run(directory, download, seed, n, partitions, device=None, class_counts=(2, 3, 4)):
    rng = np.random.default_rng(seed)
    n_fit, n_cal = 8192, 2048
    betas = np.geomspace(2, 2048, 16)
    report = dict(dataset=REPO, revision=REVISION, seed=seed, source_states=n, slopes=len(betas),
                  agreement_tolerance_bits=0.05, models={})
    reference, selection = {}, {}
    for model in (RESNET, MOBILE):
        arrays, labels = {}, {}
        for split, count in (("train", n_fit + n_cal), ("test", n)):
            shard = shard_file(directory, model, split, download)
            meta = metadata(shard, model)
            if split not in reference:
                reference[split] = meta
                selection[split] = rng.choice(len(meta), count, replace=False)
            chosen = meta[selection[split]]
            arrays[split], labels[split] = features(shard, chosen[:, 0]), chosen[:, 2]
        train, test = arrays["train"], arrays["test"]
        fit, cal = slice(0, n_fit), slice(n_fit, None)
        tasks = []
        for classes in class_counts:
            for part in range(partitions):
                groups = np.random.default_rng(1000 * classes + part).permutation(20) % classes
                y_tr, y_te = groups[labels["train"]], groups[labels["test"]]
                head, temp = calibrated_head(train[fit], y_tr[fit], train[cal], y_tr[cal], classes)
                q = np.maximum(head(test), 1e-9)
                q /= q.sum(1, keepdims=True)
                full, half = curve(q, betas, device), curve(q[: n // 2], betas, device)
                tasks.append(dict(classes=classes, rank_bound=classes - 1, partition=groups.tolist(),
                                  temperature=temp, accuracy=float(np.mean(q.argmax(1) == y_te)),
                                  predicted_bits_per_halving=(classes - 1) / 2,
                                  measured=slope_summary(full, half, 0.05), curve=full, curve_half=half))
                m = tasks[-1]["measured"]
                print(model[:9], classes, part, round(tasks[-1]["accuracy"], 3), m and
                      {k: np.round(v, 3).tolist() for k, v in m.items()}, flush=True)
        report["models"][model] = tasks
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--source", type=int, default=4096)
    parser.add_argument("--partitions", type=int, default=3)
    parser.add_argument("--device", help="run the Blahut iteration with PyTorch on this device, e.g. cuda")
    parser.add_argument("--classes", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "rank_cifar100.json")
    args = parser.parse_args()
    with threadpool_limits(limits=4):
        result = run(args.data, args.download, args.seed, args.source, args.partitions, args.device, tuple(args.classes))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n", encoding="utf-8")
