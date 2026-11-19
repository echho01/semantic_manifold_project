"""Does a learned coding surface remove the isometry defect of raw features?

The evaluation cost is the joint cost of the letter in its lifted form,
C^g_ij = || z_i - z_j ||^2 with z = (2 sqrt(q(u)), sqrt(eta) u), for a frozen calibrated head q.
A coder acts on a surface F through the Gibbs kernel W_ij ~ r_j exp(-s || F(u_i) - F(a_j) ||^2)
with one global scale s and masses r that are Blahut-optimal for the surface's own cost. Its
information is compared, at equal evaluation distortion D_g, with the Blahut optimum of C^g.

  raw      F(u) = u. The surface carries the feature geometry and is mismatched to C^g.
  learned  F(u) = (c u, MLP(u)), trained on training images to minimize the excess Lagrangian
           sum_s [ I(W_s) + s D_g(W_s) ] over a grid of scales. It sees features only.
  lift     F(u) = z, the isometric surface. Its excess is zero by construction.

All curves are evaluated on held-out official test images with anchors from the training set.
Rates are bits per image. Operating points are achievable values of finite problems.

    python experiments/raw_surface.py --data data/semasia --device cuda --output results/raw_surface.json
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp, softmax
from sklearn.linear_model import LogisticRegression
from torch import nn

from semasia_io import MOBILE, REPO, RESNET, REVISION, features, metadata, shard_file

LN2 = math.log(2)


def head(x_fit, y_fit, x_cal, y_cal):
    clf = LogisticRegression(C=10.0, max_iter=400, tol=1e-6).fit(x_fit, y_fit)
    logits = clf.decision_function(x_cal)
    nll = lambda t: float(np.mean(logsumexp(logits / math.exp(t), axis=1)
                                  - (logits / math.exp(t))[np.arange(len(y_cal)), y_cal]))
    temp = math.exp(minimize_scalar(nll, bounds=(math.log(0.1), math.log(10)), method="bounded").x)
    return (lambda x: softmax(clf.decision_function(x) / temp, axis=1)), temp


def squared(x, y):
    return ((x * x).sum(1)[:, None] + (y * y).sum(1)[None, :] - 2 * x @ y.T).clamp_min(0)


def blahut(cost, s, iters=400):
    """Masses that are optimal for `cost` at scale s (nats), uniform source."""
    kernel = torch.exp(-s * (cost - cost.min(1, keepdim=True).values))
    r = torch.full((cost.shape[1],), 1.0 / cost.shape[1], dtype=cost.dtype, device=cost.device)
    for _ in range(iters):
        r = r * (kernel / (kernel @ r)[:, None]).mean(0)
        r = r / r.sum()
    w = kernel * r
    return w / w.sum(1, keepdim=True)


def point(w, evaluation):
    marginal = w.mean(0)
    info = (w * (w.clamp_min(1e-300).log() - marginal.clamp_min(1e-300).log())).sum(1).mean() / LN2
    return float(info), float((w * evaluation).sum(1).mean())


def curve(coding, evaluation, scales):
    pts = [point(blahut(coding, float(s)), evaluation) for s in scales]
    return sorted(pts, key=lambda p: p[1])          # increasing distortion


def rate_at(pts, d):
    """Least information among the computed points whose distortion does not exceed d."""
    ok = [p[0] for p in pts if p[1] <= d]
    return min(ok) if ok else None


class Surface(nn.Module):
    def __init__(self, dim, width, out):
        super().__init__()
        self.scale = nn.Parameter(torch.tensor(1.0))
        self.net = nn.Sequential(nn.Linear(dim, width), nn.GELU(), nn.Linear(width, out))

    def forward(self, u):
        return torch.cat((self.scale * u, self.net(u)), 1)


def train_surface(u, anchors, evaluation, scales, seed, steps, width, out):
    torch.manual_seed(seed)
    net = Surface(u.shape[1], width, out).to(u.device).double()
    log_mass = nn.Parameter(torch.zeros(len(scales), len(anchors), dtype=u.dtype, device=u.device))
    opt = torch.optim.Adam([*net.parameters(), log_mass], lr=2e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
    s = torch.as_tensor(scales, dtype=u.dtype, device=u.device)
    for _ in range(steps):
        coding = squared(net(u), net(anchors))
        prior = log_mass.log_softmax(1)                                   # [S, M]
        logw = (prior[:, None, :] - s[:, None, None] * coding[None]).log_softmax(2)
        w = logw.exp()
        # I(W) <= E KL(W_i || r); adding s * D_g gives the Lagrangian of the evaluation cost.
        loss = ((w * (logw - prior[:, None, :])).sum(2) + s[:, None] * (w * evaluation[None]).sum(2)).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()
    return net.eval()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--eta", type=float, default=0.1, help="weight of u in the lift; the joint metric has half of it")
    ap.add_argument("--model", choices=("resnet", "mobilenet"), default="resnet")
    ap.add_argument("--anchors", type=int, default=256)
    ap.add_argument("--test", type=int, default=4096)
    ap.add_argument("--steps", type=int, default=3000)
    ap.add_argument("--width", type=int, default=256)
    ap.add_argument("--out-dim", type=int, default=32)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1115, 1116, 1117])
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    dev = torch.device(a.device)
    model = RESNET if a.model == "resnet" else MOBILE
    rng = np.random.default_rng(20261118)
    n_fit, n_cal = 8192, 2048
    arrays, labels = {}, {}
    for split, count in (("train", n_fit + n_cal), ("test", a.test)):
        shard = shard_file(a.data, model, split, False)
        meta = metadata(shard, model)
        chosen = meta[rng.choice(len(meta), count, replace=False)]
        arrays[split], labels[split] = features(shard, chosen[:, 0]), chosen[:, 2]
    train, test = arrays["train"], arrays["test"]
    q, temp = head(train[:n_fit], labels["train"][:n_fit], train[n_fit:], labels["train"][n_fit:])
    t = lambda x: torch.as_tensor(x, dtype=torch.float64, device=dev)
    lift = lambda x: torch.cat((2 * t(q(x)).sqrt(), math.sqrt(a.eta) * t(x)), 1)
    u_fit, u_test, u_anchor = t(train[:n_fit]), t(test), t(train[: a.anchors])
    z_fit, z_test, z_anchor = lift(train[:n_fit]), lift(test), lift(train[: a.anchors])
    eval_fit, eval_test = squared(z_fit, z_anchor), squared(z_test, z_anchor)

    reference = curve(eval_test, eval_test, np.geomspace(0.3, 300, 120))
    raw_scales = np.geomspace(0.3, 3000, 120)
    surfaces = {"raw": [curve(squared(u_test, u_anchor), eval_test, raw_scales)]}
    train_scales = np.geomspace(1, 200, 10)
    surfaces["learned"] = []
    with torch.enable_grad():
        for seed in a.seeds:
            net = train_surface(u_fit, u_anchor, eval_fit, train_scales, seed, a.steps, a.width, a.out_dim)
            with torch.no_grad():
                surfaces["learned"].append(curve(squared(net(u_test), net(u_anchor)), eval_test, raw_scales))
    floors = {name: [min(p[1] for p in c) for c in curves] for name, curves in surfaces.items()}
    floors["reference"] = [min(p[1] for p in reference)]
    print(json.dumps(floors))
    # Targets inside the range that the raw surface can reach, and below it.
    targets = [1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5]
    rows = []
    for d in targets:
        ref = rate_at(reference, d)
        row = dict(distortion=d, reference_bits=ref)
        for name, curves in surfaces.items():
            ex = [rate_at(c, d) for c in curves]
            row[name + "_excess_bits"] = [None if e is None else e - ref for e in ex]
        rows.append(row)
        print(json.dumps(row))
    out = dict(dataset=REPO, revision=REVISION, model=model, lift_weight=a.eta, eta=a.eta / 2, anchors=a.anchors, test=a.test,
               fit=n_fit, temperature=temp, seeds=a.seeds, steps=a.steps, train_scales=train_scales.tolist(),
               head_accuracy=float(np.mean(q(test).argmax(1) == labels["test"])),
               least_distortion=floors, rows=rows, reference_curve=reference, curves=surfaces)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    with torch.no_grad():
        main()
