"""Measure the component of a trained coding cost that cannot be absorbed into anchor masses."""
import argparse
import hashlib
import json
from pathlib import Path

import torch

from link import information, post, test_source
from semasia_curve import Surface, phi, squared


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--checkpoints", type=Path, nargs="+", required=True)
    parser.add_argument("--count", type=int, default=2048)
    parser.add_argument("--seed", type=int, default=31415)
    parser.add_argument("--split", choices=("validation", "test"), default="test")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    source_hash, reference = None, None
    for path in args.checkpoints:
        ck = torch.load(path, map_location="cpu", weights_only=True)
        if ck.get("kind") != "curve":
            raise ValueError("expected a coding-surface checkpoint")
        x, q, _, ids = test_source(args, ck, torch.device("cpu"))
        source_ids = hashlib.sha256(ids.cpu().numpy().tobytes()).hexdigest()
        bx = ck["anchors"]
        u, b = phi(x, q, ck["eta"]), phi(bx, post(bx, ck), ck["eta"])
        net = Surface(**ck["architecture"])
        net.load_state_dict(ck["state"])
        with torch.no_grad():
            base = squared(u, b)
            if reference is not None and (source_ids != source_hash or
                                          not torch.equal(base, reference)):
                raise ValueError("the paired experiment needs identical evaluation points and costs")
            reference, source_hash = base, source_ids
            learned = squared(net(u), net(b))
            delta = ck["scale"] * (learned.double() - base.double())
            residual = delta - delta.mean(0) - delta.mean(1, keepdim=True) + delta.mean()
            w = net.kernel(u, b, ck["scale"])
            distortion = float((w * base).sum(1).mean())
            rate = float(information(w))
            records.append(dict(mode=net.mode, checkpoint=path.name,
                                checkpoint_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                rmse_interaction_logit=float(residual.square().mean().sqrt()),
                                information_bits=rate, fixed_distortion=distortion,
                                rate_plus_lambda_distortion=rate + ck["lam"] * distortion,
                                scale=ck["scale"], lam=ck["lam"], epochs=ck["epochs"],
                                anchors=ck["book"], fit=ck["fit"], cal=ck["cal"]))
    result = dict(split=args.split, count=len(ids), seed=args.seed, id_sha256=source_hash,
                  revision=ck["revision"], model=ck["model"], task=ck["task"], modes=records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", "utf-8")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
