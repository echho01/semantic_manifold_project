"""Compute a resource converse for ONE fixed dictionary shared by all test images (complex AWGN)."""

import argparse, json, math
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import DataLoader
from images import Images
from structure import merge_values
from topology_jscc import Codec


def h(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ("data", "checkpoint", "output"):
        ap.add_argument("--" + k, type=Path, required=True)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--split", choices=("val", "test"), default="test")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--seed", type=int, default=31415)
    ap.add_argument("--queries", type=int, default=4)
    ap.add_argument("--snr", type=float, default=0)
    ap.add_argument(
        "--dictionary",
        type=Path,
        help="predeclared JSON array of [flat vertex a, flat vertex b, threshold]",
    )
    ap.add_argument("--alpha", type=float, default=0.05)
    a = ap.parse_args()
    dev = torch.device(a.device)
    if not 0 < a.alpha < 1:
        raise ValueError("0 < alpha < 1 required")
    ck = torch.load(a.checkpoint, map_location=dev, weights_only=True)
    cfg = ck["cfg"]
    s = cfg["size"]
    ds = Images(a.data, cfg["dataset"], a.split, s, a.limit, a.seed)
    rng = np.random.default_rng(a.seed)
    if a.dictionary:
        queries = np.asarray(json.loads(a.dictionary.read_text("utf8")), float)
    else:
        if not 1 <= a.queries <= 8:
            raise ValueError("choose 1 to 8 fixed queries")
        queries = np.column_stack(
            (rng.integers(s * s, size=(a.queries, 2)), rng.uniform(0.2, 0.8, a.queries))
        )
    if queries.ndim != 2 or queries.shape[1] != 3 or not 1 <= len(queries) <= 8:
        raise ValueError("1 to 8 triples required")
    pairs = queries[:, :2].astype(int)
    t = queries[:, 2]
    if (
        (queries[:, :2] != pairs).any()
        or pairs.min() < 0
        or pairs.max() >= s * s
        or not np.all((0 <= t) & (t <= 1))
    ):
        raise ValueError("invalid vertices/thresholds")
    net = Codec(s, cfg["rho"], cfg["hidden"], cfg["groups"], cfg["power"]).to(dev)
    net.load_state_dict(ck["state"])
    net.eval()
    torch.manual_seed(a.seed + 7)
    truth = []
    received = []
    with torch.no_grad():
        for x, _ in DataLoader(ds, batch_size=32):
            x = x.to(dev)
            y, _, _ = net(x, "awgn", a.snr)
            for ref, r in zip(x[:, 0], y[:, 0]):
                truth.append(merge_values(1 - ref, pairs).cpu().numpy() <= t)
                received.append(merge_values(1 - r, pairs).cpu().numpy() <= t)
    z = np.array(truth)
    errors = (z != received).mean(0)
    n, d = z.shape
    K = 2**d
    pr = np.bincount(z @ (1 << np.arange(d)), minlength=K) / n
    H = float(-(pr[pr > 0] * np.log2(pr[pr > 0])).sum())
    # TV=sup_A |P(A)-Pn(A)|; a union bound over subsets and Hoeffding (allocate alpha/2).
    tv = math.sqrt((K * math.log(2) + math.log(4 / a.alpha)) / (2 * n))
    Hlo = 0 if tv >= 1 - 1 / K else max(0, H - float(h(tv)) - tv * math.log2(K - 1))
    e_hi = np.minimum(0.5, errors + math.sqrt(math.log(2 * d / a.alpha) / (2 * n)))
    lower = max(0.0, Hlo - float(h(e_hi).sum()))
    out = dict(
        images=n,
        split=a.split,
        ids_sha256=ds.sha,
        dictionary=queries.tolist(),
        snr_db=a.snr,
        H_plugin_bits=H,
        error_per_query=errors.tolist(),
        plugin_converse_bits=max(0.0, H - float(h(errors).sum())),
        lower_bits_confidence=lower,
        confidence_if_iid=1 - a.alpha,
        actual_k=net.k,
        lower_awgn_uses_if_iid=lower / math.log2(1 + 10 ** (a.snr / 10)),
        assumption="Fixed dictionary chosen independently of this split; confidence assumes independent images.",
    )
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2) + "\n", "utf8")
    print(a.output)


if __name__ == "__main__":
    main()
