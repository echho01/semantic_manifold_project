"""Paired summary of the fine-tuned arms: per-image differences to the MSE arm, averaged over seeds."""
import gzip, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

root = Path(sys.argv[1]); seeds = (1115, 1116, 1117); arms = ("mse", "query", "proposed", "topojscc")
cols = ["image", "channel", "snr_db", "repeat", "mse", "hybrid_query", "uniform_query", "w2_h0", "w2_h1"]
data = {}
for arm in arms:
    for s in seeds:
        f = root / f"{arm}_{s}_paired.jsonl.gz"
        if not f.exists():
            continue
        table = defaultdict(dict)
        with gzip.open(f, "rt") as fh:
            for line in fh:
                r = dict(zip(cols, json.loads(line)))
                table[(r["channel"], r["snr_db"])][r["image"]] = (r["mse"], r["hybrid_query"], r["uniform_query"], r["w2_h0"])
        data[arm, s] = table
out = []
conditions = sorted({c for t in data.values() for c in t})
for cond in conditions:
    row = dict(channel=cond[0], snr_db=cond[1], arms={})
    for arm in arms:
        ss = [s for s in seeds if (arm, s) in data and ("mse", s) in data]
        if not ss:
            continue
        images = sorted(data[arm, ss[0]][cond])
        a = np.array([[data[arm, s][cond][i] for i in images] for s in ss])        # [seed, image, metric]
        b = np.array([[data["mse", s][cond][i] for i in images] for s in ss])
        psnr = -10 * np.log10(a[..., 0].mean(1))
        diff = (a[..., 1] - b[..., 1]).mean(0)                                      # per image, mean over seeds
        half = 1.96 * diff.std(ddof=1) / np.sqrt(len(diff))
        row["arms"][arm] = dict(seeds=len(ss), images=len(images), psnr_db=float(psnr.mean()),
                                psnr_by_seed=psnr.round(3).tolist(),
                                hybrid_risk=float(a[..., 1].mean()), hybrid_risk_by_seed=a[..., 1].mean(1).round(5).tolist(),
                                uniform_risk=float(a[..., 2].mean()), w2_h0=float(a[..., 3].mean()),
                                risk_minus_mse=float(diff.mean()), ci95=[float(diff.mean() - half), float(diff.mean() + half)],
                                relative_change=float(diff.mean() / b[..., 1].mean()))
    out.append(row)
    print(cond, {k: (round(v["psnr_db"], 2), round(v["hybrid_risk"], 4), round(100 * v["relative_change"], 1),
                     [round(100 * x / row["arms"]["mse"]["hybrid_risk"], 1) for x in v["ci95"]]) for k, v in row["arms"].items()})
json.dump(out, open(root / "summary.json", "w"), indent=1)
