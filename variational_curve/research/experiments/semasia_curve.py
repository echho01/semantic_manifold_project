"""Bounded residual coding surface and a separate mass per fixed anchor."""

import math
import torch
from torch import nn
from link import (
    arguments,
    device,
    train_source,
    test_source,
    post,
    report,
    save,
    LN2,
    information,
    kl_table,
)


def phi(x, q, eta):
    return torch.cat((2 * q.sqrt(), math.sqrt(eta) * x), 1)


def squared(x, y):
    return ((x * x).sum(1)[:, None] + (y * y).sum(1)[None, :] - 2 * x @ y.T).clamp_min(
        0
    )


class Surface(nn.Module):
    def __init__(self, dim, rank, m, alpha, mode):
        super().__init__()
        self.a = nn.Parameter(torch.randn(rank, dim) / math.sqrt(dim))
        self.b = nn.Parameter(torch.zeros(dim, rank))
        self.log_mass = nn.Parameter(torch.zeros(m))
        self.alpha, self.mode = alpha, mode

    def forward(self, x):
        if self.mode not in ("joint", "geometry"):
            return x
        # Frobenius bounds the operator norm at every optimizer step, including on unseen features.
        A = self.a / self.a.norm().clamp_min(1)
        B = self.b / self.b.norm().clamp_min(1)
        return x + self.alpha * ((x @ A.T).tanh() @ B.T)

    def kernel(self, x, b, s):
        m = (
            self.log_mass
            if self.mode in ("density", "joint")
            else torch.zeros_like(self.log_mass)
        )
        return (m - s * squared(self(x), self(b))).softmax(1)


def main():
    ap, tr, _ = arguments(__doc__)
    for name, typ, val in [
        ("epochs", int, 30),
        ("batch", int, 256),
        ("rank", int, 12),
        ("alpha", float, 0.6),
        ("eta", float, 0.1),
        ("scale", float, 5.0),
        ("lam", float, 3.0),
        ("bend", float, 0.05),
        ("lr", float, 0.01),
    ]:
        tr.add_argument("--" + name, type=typ, default=val)
    tr.add_argument(
        "--mode", choices=("fixed", "density", "geometry", "joint"), default="joint"
    )
    a = ap.parse_args()
    dev = device(a.device)
    if a.command == "train":
        if not (
            0 <= a.alpha < 1
            and a.rank > 0
            and a.eta > 0
            and a.scale > 0
            and a.lam >= 0
            and a.bend >= 0
            and a.epochs > 0
            and a.batch > 0
            and a.lr > 0
        ):
            raise ValueError("invalid optimizer or surface settings")
        torch.manual_seed(a.seed)
        ck, x, q = train_source(a, dev)
        bq = post(ck["anchors"].to(dev), ck)
        u = phi(x, q, a.eta)
        b = phi(ck["anchors"].to(dev), bq, a.eta)
        net = Surface(u.shape[1], a.rank, len(b), a.alpha, a.mode).to(dev)
        # Optimized auxiliary marginal gives an upper bound on I, exact at the true output marginal.
        aux = nn.Parameter(torch.zeros(len(b), device=dev))
        opt = torch.optim.Adam([*net.parameters(), aux], lr=a.lr)
        history = []
        for epoch in range(a.epochs):
            for ix in torch.randperm(len(u), device=dev).split(a.batch):
                z = u[ix]
                c = squared(z, b)
                W = net.kernel(z, b, a.scale)
                rate = (W * (W.clamp_min(1e-30).log() - aux.log_softmax(0))).sum(
                    1
                ).mean() / LN2
                # Nearest distinct anchor; chord midpoint tests the extension of F, not a Willmore estimator.
                d = c.detach().masked_fill(c.detach() < 1e-6, float("inf"))
                v = b[d.argmin(1)]
                mid = (z + v) / 2
                bend = (net(z) + net(v) - 2 * net(mid)).square().sum(1).mean()
                loss = rate + a.lam * (W * c).sum(1).mean() + a.bend * bend
                opt.zero_grad()
                loss.backward()
                opt.step()
            with torch.no_grad():
                W = torch.cat(
                    [net.kernel(block, b, a.scale) for block in u.split(a.batch)]
                )
                history.append(
                    dict(
                        epoch=epoch + 1,
                        information_bits=float(information(W)),
                        fixed_distortion=float((W * squared(u, b)).sum(1).mean()),
                        kl_bits=float((W * kl_table(q, bq)).sum(1).mean()),
                    )
                )
        ck.update(
            kind="curve",
            architecture=dict(
                dim=u.shape[1], rank=a.rank, m=len(b), alpha=a.alpha, mode=a.mode
            ),
            state={k: v.cpu() for k, v in net.state_dict().items()},
            eta=a.eta,
            scale=a.scale,
            lam=a.lam,
            bend=a.bend,
            epochs=a.epochs,
            batch=a.batch,
            lr=a.lr,
            history=history,
        )
        print(history[-1])
        save(a, ck)
    else:
        ck = torch.load(a.checkpoint, map_location="cpu", weights_only=True)
        if ck.get("kind") != "curve":
            raise ValueError("not a surface checkpoint")
        x, q, y, ids = test_source(a, ck, dev)
        bx = ck["anchors"].to(dev)
        bq = post(bx, ck)
        u, b = phi(x, q, ck["eta"]), phi(bx, bq, ck["eta"])
        net = Surface(**ck["architecture"]).to(dev)
        net.load_state_dict(ck["state"])
        net.eval()
        with torch.no_grad():
            W = torch.cat([net.kernel(z, b, ck["scale"]) for z in u.split(512)])
            report(
                a,
                ck,
                q,
                y,
                ids,
                W,
                bq,
                dict(
                    mode=net.mode,
                    lambda_fixed=ck["lam"],
                    expected_fixed_distortion=float((W * squared(u, b)).sum(1).mean()),
                    train_history=ck["history"],
                ),
                cosine=x @ bx.T,
                fixed=squared(u, b),
            )


if __name__ == "__main__":
    main()
