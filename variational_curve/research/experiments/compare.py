"""Paired comparison at a fixed physical operating point; resample independent test groups."""

import argparse, gzip, json
from pathlib import Path
import numpy as np


def observations(path, channel, snr, metric):
    j = json.loads(path.read_text("utf8"))
    if "k_complex" in j:
        cols = j["paired_columns"]
        ix = cols.index(metric)
        x = {}
        groups = {}
        with gzip.open(path.parent / j["paired_file"], "rt", encoding="utf8") as handle:
            for line in handle:
                row = json.loads(line)
                key, ch, db, rep = row[:4]
                if ch != channel or db != snr:
                    continue
                x.setdefault(key, []).append(float(row[ix]))
                groups[key] = (
                    key.split("/")[1] if j["config"]["dataset"] == "omniglot" else key
                )
        context = [
            j["config"]["dataset"],
            j["split"],
            j["config"]["size"],
            j["actual_cbr"],
            j["id_sha256"],
            j["queries"],
            j["seed"],
        ]
    else:
        if metric not in ("kl_bits", "accuracy"):
            raise ValueError("use kl_bits or accuracy for SEMASIA paired arrays")
        code = channel if channel == "clean" else channel + "_" + str(float(snr))
        p = path.parent / (path.stem + "_paired") / (code + ".npz")
        with np.load(p, allow_pickle=False) as d:
            arr = d["kl_bits_accuracy"][:, :, 0 if metric == "kl_bits" else 1]
            x = {str(i): arr[:, col].tolist() for col, i in enumerate(d["id"])}
        groups = {k: k for k in x}
        context = [
            j["dataset"],
            j["model"],
            j["task"],
            j["test_split"],
            j["uses_per_sample"],
            j["test_id_sha256"],
            j["test_seed"],
        ]
    if not x:
        raise ValueError("specified condition has no paired observations")
    return context, x, groups


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("first", type=Path)
    p.add_argument("second", type=Path)
    p.add_argument("--channel", required=True)
    p.add_argument("--snr", type=float)
    p.add_argument("--metric", required=True)
    p.add_argument("--draws", type=int, default=2000)
    p.add_argument("--seed", type=int, default=31415)
    a = p.parse_args()
    c, x, groups = observations(a.first, a.channel, a.snr, a.metric)
    k, y, _ = observations(a.second, a.channel, a.snr, a.metric)
    if c != k or x.keys() != y.keys():
        raise ValueError(
            "paired sources, group identities, evaluation or physical budgets differ"
        )
    if any(len(x[i]) != len(y[i]) for i in x):
        raise ValueError("unequal channel repeats")
    if a.draws < 2:
        raise ValueError("need at least two bootstrap draws")
    d = {}
    for i in x:
        d.setdefault(groups[i], []).append(np.mean(np.array(x[i]) - np.array(y[i])))
    sums = np.array([sum(vals) for vals in d.values()])
    counts = np.array([len(vals) for vals in d.values()])
    rng = np.random.default_rng(a.seed)
    boot = []
    for _ in range(a.draws):
        ix = rng.integers(len(d), size=len(d))
        boot.append(sums[ix].sum() / counts[ix].sum())
    print(
        json.dumps(
            dict(
                metric=a.metric,
                first_minus_second=sums.sum() / counts.sum(),
                groups=len(d),
                images=len(x),
                percentile95=np.quantile(boot, [0.025, 0.975]).tolist(),
                scope="test groups conditional on these two fitted checkpoints",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
