"""Common experimental contract, vendored locally so each paper travels independently."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from sklearn.linear_model import LogisticRegression
from threadpoolctl import threadpool_limits

from semasia_io import MOBILE, RESNET, REVISION, REPO, shard_file, metadata, features

LN2 = math.log(2)


def arguments(description):
    ap = argparse.ArgumentParser(description=description)
    sp = ap.add_subparsers(dest="command", required=True)
    tr, te = sp.add_parser("train"), sp.add_parser("test")
    for sub in (tr, te):
        sub.add_argument("--data", type=Path, required=True)
        sub.add_argument(
            "--download",
            action="store_true",
            help="fetch and SHA-256 check missing official shards",
        )
        sub.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
        sub.add_argument("--checkpoint", type=Path, required=True)
    tr.add_argument("--model", choices=(RESNET, MOBILE), default=RESNET)
    tr.add_argument("--task", choices=("fine", "coarse"), default="coarse")
    tr.add_argument("--fit", type=int, default=8192)
    tr.add_argument("--cal", type=int, default=2048)
    tr.add_argument(
        "--book",
        type=int,
        default=64,
        help="power-of-two alphabet; first training IDs are anchors",
    )
    tr.add_argument("--seed", type=int, default=1115)
    te.add_argument(
        "--split",
        choices=("test", "validation"),
        default="test",
        help="validation is unused official-train IDs",
    )
    te.add_argument("--count", type=int, default=2048)
    te.add_argument("--seed", type=int, default=31415, help="test ID and channel seed")
    te.add_argument("--ebn0", type=float, nargs="+", default=[0.0, 5.0, 10.0, 15.0])
    te.add_argument(
        "--channels",
        nargs="+",
        choices=("clean", "awgn", "rayleigh"),
        default=["clean", "awgn", "rayleigh"],
    )
    te.add_argument("--repeats", type=int, default=5)
    te.add_argument("--output", type=Path, required=True)
    return ap, tr, te


def device(name):
    return torch.device(
        ("cuda" if torch.cuda.is_available() else "cpu") if name == "auto" else name
    )


def draw(args, split, model, task, n, seed, exclude=()):
    shard = shard_file(args.data, model, split, args.download)
    meta = metadata(shard, model)
    if len(exclude):
        meta = meta[~np.isin(meta[:, 0], np.asarray(exclude))]
    if not 0 < n <= len(meta):
        raise ValueError(f"sample count must be in [1,{len(meta)}]")
    selected = meta[np.random.default_rng(seed).permutation(len(meta))[:n]]
    return (
        torch.tensor(features(shard, selected[:, 0]), dtype=torch.float32),
        torch.tensor(selected[:, 1 if task == "fine" else 2], dtype=torch.long),
        torch.tensor(selected[:, 0].astype(np.int64)),
    )


def post(x, ck):
    logits = (x @ ck["weight"].to(x) + ck["bias"].to(x)) / ck["temperature"]
    q = logits.softmax(1).clamp_min(1e-12)
    return q / q.sum(1, keepdim=True)


def train_source(a, dev):
    if a.book < 1 or a.book & (a.book - 1) or a.book >= a.fit or a.cal < 1:
        raise ValueError(
            "book must be a power of two smaller than fit; calibration must be nonempty"
        )
    x, y, ids = draw(a, "train", a.model, a.task, a.fit + a.cal, a.seed)
    with threadpool_limits(limits=1):
        h = LogisticRegression(C=10, max_iter=600).fit(
            x[: a.fit].numpy(), y[: a.fit].numpy()
        )
    classes = 20 if a.task == "coarse" else 100
    if not np.array_equal(h.classes_, np.arange(classes)):
        raise ValueError("increase --fit to cover every task class")
    z = h.decision_function(x[a.fit :].numpy())
    target = y[a.fit :].numpy()

    def loss(b):
        u = z / np.exp(b)
        return (logsumexp(u, axis=1) - u[np.arange(len(u)), target]).mean()

    opt = minimize_scalar(loss, bounds=(np.log(0.1), np.log(10.0)), method="bounded")
    if not opt.success:
        raise RuntimeError("temperature calibration failed")
    ck = dict(
        revision=REVISION,
        dataset=REPO,
        model=a.model,
        task=a.task,
        seed=a.seed,
        fit=a.fit,
        cal=a.cal,
        weight=torch.tensor(h.coef_.T, dtype=torch.float32),
        bias=torch.tensor(h.intercept_, dtype=torch.float32),
        temperature=float(np.exp(opt.x)),
        fit_ids=ids[: a.fit],
        cal_ids=ids[a.fit :],
        anchor_ids=ids[: a.book],
        anchors=x[: a.book].clone(),
        book=a.book,
    )
    u = x[: a.fit].to(dev)
    q = post(u, ck)
    return ck, u, q


def test_source(a, ck, dev):
    if (ck["revision"], ck["dataset"]) != (REVISION, REPO):
        raise ValueError("checkpoint/data revision mismatch")
    split = "test" if a.split == "test" else "train"
    exclude = (
        () if split == "test" else torch.cat((ck["fit_ids"], ck["cal_ids"])).numpy()
    )
    x, y, ids = draw(a, split, ck["model"], ck["task"], a.count, a.seed, exclude)
    x = x.to(dev)
    return x, post(x, ck), y.to(dev), ids


def kl_table(q, b):
    return ((q * q.log()).sum(1, keepdim=True) - q @ b.log().T).clamp_min(0) / LN2


def information(W):
    mu = W.mean(0)
    return (W * (W.clamp_min(1e-30).log() - mu.clamp_min(1e-30).log())).sum(
        1
    ).mean() / LN2


def transmit(j, width, ch, snr, gen):
    if width == 0:
        return j, 0
    shifts = torch.arange(width - 1, -1, -1, device=j.device)
    bits = (j[:, None] >> shifts) & 1
    if ch == "clean":
        return (bits << shifts).sum(1), 0
    symbols = 2 * bits.float() - 1
    # Unit-energy BPSK; Rayleigh block amplitude is common to one index, E[a^2]=1.
    gain = (
        1
        if ch == "awgn"
        else torch.sqrt(
            -torch.rand((len(j), 1), device=j.device, generator=gen)
            .clamp_min(1e-12)
            .log()
        )
    )
    noise = torch.randn(symbols.shape, device=j.device, generator=gen) * math.sqrt(
        0.5 * 10 ** (-snr / 10)
    )
    detected = (gain * symbols + noise >= 0).long()
    return (detected << shifts).sum(1), int((bits != detected).sum())


@torch.no_grad()
def report(a, ck, q, y, ids, W, b, extra=None, cosine=None, fixed=None):
    if a.repeats < 1:
        raise ValueError("--repeats must be positive")
    width = (len(b) - 1).bit_length()
    bits = kl_table(q, b)
    gen = torch.Generator(device=q.device)
    output = dict(
        dataset=REPO,
        revision=REVISION,
        model=ck["model"],
        task=ck["task"],
        test_split=("official test" if a.split == "test" else "unused official train"),
        test_id_sha256=hashlib.sha256(ids.numpy().astype("<i8").tobytes()).hexdigest(),
        count=len(q),
        checkpoint=str(a.checkpoint),
        test_seed=a.seed,
        book=len(b),
        uses_per_sample=width,
        information_bits=float(information(W)),
        expected_clean_kl_bits=float((W * bits).sum(1).mean()),
        controls=[],
        head_accuracy=float((q.argmax(1) == y).float().mean()),
        head_nll_bits=float(-q[torch.arange(len(y), device=y.device), y].log2().mean()),
        links=[],
        extra=extra or {},
    )
    for ch in a.channels:
        for snr in [None] if ch == "clean" else a.ebn0:
            rows = []
            control_rows = {}
            errors = []
            for rep in range(a.repeats):
                gen.manual_seed(a.seed + 1009 * rep)
                j = torch.multinomial(W, 1, generator=gen).squeeze(1)
                gen.manual_seed(a.seed + 123456 + 1009 * rep)
                recv, wrong = transmit(j, width, ch, snr, gen)
                p = b[recv]
                rows.append(
                    [
                        float(bits.gather(1, recv[:, None]).mean()),
                        float((p.argmax(1) == y).float().mean()),
                        float(
                            -p[torch.arange(len(y), device=y.device), y].log2().mean()
                        ),
                        float((j != recv).float().mean()),
                        wrong / max(1, len(y) * width),
                    ]
                )
                if fixed is not None:
                    rows[-1].append(float(fixed.gather(1, recv[:, None]).mean()))
                errors.append(
                    torch.stack(
                        (
                            bits.gather(1, recv[:, None]).squeeze(1),
                            (p.argmax(1) == y).float(),
                        ),
                        1,
                    )
                    .cpu()
                    .numpy()
                )
                # Deterministic operational baselines use the same bit noise for each repeat.
                anchors = len(ck["anchors"])
                j_kl = bits[:, :anchors].argmin(1)
                controls = [("nearest_KL_anchor", j_kl, anchors)]
                if cosine is not None:
                    controls.append(
                        ("nearest_cosine_anchor", cosine.argmax(1), anchors)
                    )
                controls.append(("direct_task_decision", q.argmax(1), q.shape[1]))
                for name, index, M in controls:
                    bw = (M - 1).bit_length()
                    gen.manual_seed(a.seed + 123456 + 1009 * rep)
                    got, nerr = transmit(index, bw, ch, snr, gen)
                    invalid = got >= M
                    if name == "direct_task_decision":
                        # Invalid bit patterns are erasures; conservative strict-error scoring.
                        correct = got == y
                        kl = None
                    else:
                        correct = b[got].argmax(1) == y
                        kl = float(bits.gather(1, got[:, None]).mean())
                    control_rows.setdefault(name, []).append(
                        dict(
                            uses=bw,
                            accuracy=float(correct.float().mean()),
                            index_bler=float((got != index).float().mean()),
                            erasure=float(invalid.float().mean()),
                            ber=nerr / max(1, len(index) * bw),
                            posterior_kl_bits=kl,
                        )
                    )

            names = ["posterior_kl_bits", "accuracy", "nll_bits", "index_bler", "ber"]
            if fixed is not None:
                names += ["fixed_distortion"]
            for name, runs in control_rows.items():
                output["controls"].append(
                    dict(method=name, channel=ch, ebn0_db=snr, replicates=runs)
                )
            pair_dir = a.output.parent / (a.output.stem + "_paired")
            pair_dir.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                pair_dir / (ch + ("_" + str(snr) if snr is not None else "") + ".npz"),
                id=ids.numpy(),
                kl_bits_accuracy=np.asarray(errors, dtype=np.float32),
            )
            output["links"].append(
                dict(
                    channel=ch,
                    ebn0_db=snr,
                    repeats=a.repeats,
                    **dict(zip(names, np.mean(rows, axis=0).tolist())),
                    replicates=rows,
                    columns=names,
                )
            )
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf8"
    )
    print(
        a.output,
        "I =",
        round(output["information_bits"], 5),
        "clean posterior KL (bit) =",
        round(output["expected_clean_kl_bits"], 5),
    )
    return output


def save(a, ck):
    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(ck, a.checkpoint)
    print("checkpoint:", a.checkpoint)
