"""Constrained complex JSCC and aligned connectivity learning on ordinary images or road masks."""

import argparse, gzip, hashlib, json, math
from pathlib import Path
import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader
from images import Images
from structure import batch_ph, connection, image_ph, mst_signature


class Codec(nn.Module):
    def __init__(self, size, rho, hidden, groups=1, power="uniform"):
        super().__init__()
        if size % 4 or size < 8 or not 0 < rho <= 1 or hidden < 1:
            raise ValueError("size must be >=8 and divisible by 4; 0<rho<=1; hidden>=1")
        self.k = max(1, round(size * size * rho))
        self.side = size // 4
        self.hidden = hidden
        if groups < 1 or groups > self.k or (power == "learned" and groups < 2):
            raise ValueError(
                "power groups must be in [1,k]; learning requires at least two"
            )
        self.groups = groups
        self.log_power = nn.Parameter(
            torch.zeros(groups), requires_grad=power == "learned"
        )
        channels = [1] + [hidden] * 5
        enc = []
        for i, s in enumerate([2, 2, 1, 1, 1]):
            enc.extend([nn.Conv2d(channels[i], channels[i + 1], 5, s, 2), nn.PReLU()])
        self.enc = nn.Sequential(*enc)
        self.pack = nn.Linear(hidden * self.side**2, 2 * self.k)
        self.unpack = nn.Linear(2 * self.k, hidden * self.side**2)
        dec = []
        for i, s in enumerate([1, 1, 1, 2, 2]):
            dec.extend(
                [
                    nn.ConvTranspose2d(
                        hidden, 1 if i == 4 else hidden, 5, s, 2, output_padding=s - 1
                    ),
                    nn.Sigmoid() if i == 4 else nn.PReLU(),
                ]
            )
        self.dec = nn.Sequential(*dec)

    def forward(self, x, channel, snr):
        z = self.pack(self.enc(x).flatten(1)).reshape(len(x), self.k, 2)
        chunks = z.tensor_split(self.groups, dim=1)
        lens = z.new_tensor([v.shape[1] for v in chunks])
        p = (self.log_power + lens.log()).softmax(0) * self.k / lens
        z = torch.cat(
            [
                v
                * (
                    pi / v.square().sum(2).mean(1, keepdim=True).clamp_min(1e-12)
                ).sqrt()[:, :, None]
                for v, pi in zip(chunks, p)
            ],
            1,
        )
        s = 10 ** (-snr / 10)
        n = torch.randn_like(z) * math.sqrt(s / 2)
        if channel == "awgn":
            rx = z + n
        elif channel == "rayleigh":
            h = torch.randn((len(x), 1, 2), device=x.device) / math.sqrt(2)
            a, b = h[..., 0], h[..., 1]
            re = a * z[..., 0] - b * z[..., 1] + n[..., 0]
            im = b * z[..., 0] + a * z[..., 1] + n[..., 1]
            den = a * a + b * b + s
            rx = torch.stack(((a * re + b * im) / den, (a * im - b * re) / den), 2)
        else:
            raise ValueError(channel)
        y = self.dec(
            self.unpack(rx.flatten(1)).reshape(-1, self.hidden, self.side, self.side)
        )
        return y, z.flatten(1), rx.flatten(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sp = parser.add_subparsers(dest="action", required=True)
    train, test = sp.add_parser("train"), sp.add_parser("test")
    for p in (train, test):
        p.add_argument("--data", type=Path, required=True)
        p.add_argument("--checkpoint", type=Path, required=True)
        p.add_argument("--device", choices=("cuda", "cpu", "auto"), default="auto")
        p.add_argument("--limit", type=int)
        p.add_argument("--batch", type=int, default=32)
        p.add_argument("--workers", type=int, default=0)
        p.add_argument("--seed", type=int, default=1115)
        p.add_argument("--queries", type=int, default=32)
        p.add_argument(
            "--channels",
            nargs="+",
            choices=("awgn", "rayleigh"),
            default=["awgn", "rayleigh"],
        )
        p.add_argument(
            "--snr", type=float, nargs="+", default=[0.0, 5.0, 10.0, 15.0, 20.0]
        )
    train.add_argument("--dataset", choices=("omniglot", "roads"), default="omniglot")
    train.add_argument("--groups", type=int, default=1)
    train.add_argument("--power", choices=("uniform", "learned"), default="uniform")
    train.add_argument(
        "--mode",
        choices=("mse", "topojscc", "topoae", "query", "proposed"),
        default="proposed",
    )
    for key, typ, default in [
        ("size", int, 32),
        ("rho", float, 0.4),
        ("hidden", int, 24),
        ("epochs", int, 30),
        ("lr", float, 1e-4),
        ("lambda-query", float, 0.2),
        ("lambda-mst", float, 0.02),
        ("lambda-img", float, 1e-4),
        ("lambda-lat", float, 1e-5),
        ("ph-batch", int, 12),
        ("structure-images", int, 8),
        ("warmup", float, 2.0),
    ]:
        train.add_argument("--" + key, type=typ, default=default)
    test.add_argument("--split", choices=("val", "test"), default="test")
    test.add_argument("--repeats", type=int, default=3)
    test.add_argument("--output", type=Path, required=True)
    a = parser.parse_args()
    if a.batch < 1 or a.queries < 1 or a.workers < 0:
        raise ValueError("invalid batch/query/worker count")
    dev = torch.device(
        ("cuda" if torch.cuda.is_available() else "cpu")
        if a.device == "auto"
        else a.device
    )
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    query_rng = np.random.default_rng(a.seed + 97)
    if a.action == "train":
        if (
            min(a.epochs, a.ph_batch, a.structure_images, a.lr) <= 0
            or min(a.lambda_mst, a.lambda_query, a.lambda_img, a.lambda_lat, a.warmup)
            < 0
        ):
            raise ValueError("invalid training settings")
        net = Codec(a.size, a.rho, a.hidden, a.groups, a.power).to(dev)
        ds = Images(a.data, a.dataset, "train", a.size, a.limit, a.seed)
        dl = DataLoader(ds, batch_size=a.batch, shuffle=True, num_workers=a.workers)
        opt = torch.optim.Adam(net.parameters(), lr=a.lr)
        history = []
        for ep in range(a.epochs):
            metrics = []
            for x, _ in dl:
                x = x.to(dev)
                ch = str(rng.choice(a.channels))
                snr = float(rng.choice(a.snr))
                y, z, r = net(x, ch, snr)
                mse = (y - x).square().mean()
                loss = mse
                w = 1 if a.warmup == 0 else 1 - math.exp(-(ep + 1) / a.warmup)
                sl = slice(0, a.structure_images)
                lat = slice(0, a.ph_batch)
                if a.mode == "topojscc":
                    cost = image_ph(x[sl], y[sl]).sum()
                    if len(x) > 1:
                        cost = a.lambda_img * cost + a.lambda_lat * batch_ph(
                            z[lat], r[lat]
                        )
                    else:
                        cost = a.lambda_img * cost
                    loss = loss + w * cost
                if a.mode in ("query", "proposed"):
                    loss = loss + w * a.lambda_query * connection(
                        x[sl], y[sl], a.queries, query_rng
                    )
                if a.mode in ("topoae", "proposed"):
                    loss = loss + w * a.lambda_mst * mst_signature(x[lat], r[lat])
                opt.zero_grad()
                loss.backward()
                opt.step()
                metrics.append((float(loss.detach()), float(mse.detach()), len(x)))
            arr = np.array(metrics)
            history.append(
                dict(
                    epoch=ep + 1,
                    loss=float(np.average(arr[:, 0], weights=arr[:, 2])),
                    mse=float(np.average(arr[:, 1], weights=arr[:, 2])),
                )
            )
        a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
        cfg = {k: v for k, v in vars(a).items() if k not in ("data", "checkpoint")}
        torch.save(
            dict(
                cfg=cfg, train_id_sha256=ds.sha, history=history, state=net.state_dict()
            ),
            a.checkpoint,
        )
        print(a.checkpoint, history[-1])
        return
    if a.repeats < 1:
        raise ValueError("repeats must be positive")
    ck = torch.load(a.checkpoint, map_location=dev, weights_only=True)
    cfg = ck["cfg"]
    net = Codec(cfg["size"], cfg["rho"], cfg["hidden"], cfg["groups"], cfg["power"]).to(
        dev
    )
    net.load_state_dict(ck["state"])
    net.eval()
    ds = Images(a.data, cfg["dataset"], a.split, cfg["size"], a.limit, a.seed)
    dl = DataLoader(ds, batch_size=a.batch, num_workers=a.workers)
    rows = []
    a.output.parent.mkdir(parents=True, exist_ok=True)
    paired = a.output.with_name(a.output.stem + "_paired.jsonl.gz")
    with torch.no_grad(), gzip.open(paired, "wt", encoding="utf8") as samples:
        for ch in a.channels:
            for snr in a.snr:
                runs = []
                for rep in range(a.repeats):
                    torch.manual_seed(a.seed + 10000 + rep)
                    tot = np.zeros(5)
                    for x, paths in dl:
                        x = x.to(dev)
                        y, _, _ = net(x, ch, snr)
                        per = []
                        for i, p in enumerate(paths):
                            k = int.from_bytes(
                                hashlib.sha256((str(a.seed) + p).encode()).digest()[:8],
                                "big",
                            )
                            metrics = []
                            for kind in ("hybrid", "uniform"):
                                v = float(
                                    connection(
                                        x[i : i + 1],
                                        y[i : i + 1],
                                        a.queries,
                                        np.random.default_rng(k),
                                        kind,
                                    )
                                )
                                metrics.append(v)
                            ph = image_ph(x[i : i + 1], y[i : i + 1])
                            vals = [
                                float((x[i] - y[i]).square().mean()),
                                *metrics,
                                float(ph[0]),
                                float(ph[1]),
                            ]
                            per.append(vals)
                            samples.write(json.dumps([p, ch, snr, rep, *vals]) + "\n")
                        tot += np.asarray(per).sum(0)
                    e = tot / len(ds)
                    runs.append(
                        dict(
                            mse=e[0],
                            psnr_db=-10 * math.log10(max(e[0], 1e-16)),
                            query_risk=e[1],
                            uniform_query_risk=e[2],
                            image_w2_h0=e[3],
                            image_w2_h1=e[4],
                        )
                    )
                rows.append(
                    dict(
                        channel=ch,
                        snr_db=snr,
                        replicates=runs,
                        mean={k: float(np.mean([r[k] for r in runs])) for k in runs[0]},
                    )
                )
    lens = np.array([len(z) for z in np.array_split(np.arange(net.k), net.groups)])
    raw = net.log_power.detach().double().cpu().numpy()
    power = np.exp(raw - raw.max())
    power *= net.k / (power @ lens)
    out = dict(
        paired_file=paired.name,
        paired_columns=[
            "image",
            "channel",
            "snr_db",
            "repeat",
            "mse",
            "hybrid_query",
            "uniform_query",
            "w2_h0",
            "w2_h1",
        ],
        power_per_use=power.tolist(),
        group_sizes=lens.tolist(),
        split=a.split,
        images=len(ds),
        id_sha256=ds.sha,
        train_id_sha256=ck["train_id_sha256"],
        seed=a.seed,
        queries=a.queries,
        config=cfg,
        k_complex=net.k,
        actual_cbr=net.k / cfg["size"] ** 2,
        conditions=rows,
    )
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n", "utf8")
    print(a.output)


if __name__ == "__main__":
    main()
