"""Three controls on the allocation rules and the converse.

  analog    Does the p ~ c^(2/3) rule help at noise levels far outside the half-gap where the
            first variation is exact? Gray 16x16 glyph fields, exact all-pairs relation risk,
            common Gaussian noise for uniform and allocated power.
  shared    How much of the BPSK gain needs the powers to be computed from the image itself?
            One power map from the mean influence of training-alphabet glyphs is applied to
            every held-out glyph.
  converse  A source on which the query bound is not trivial: the eight symmetry images of a
            glyph, equiprobable, with a fixed dictionary of three connectivity queries. The
            state index is sent by BPSK and the bound is compared with the channel uses.

    python experiments/allocation_controls.py --data data/omniglot --output results/allocation_controls.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

from orbit_rdf import orbit
from relation_variation import (bpsk_allocation, components, critical_loads, disagreement, influence_map,
                                load_glyphs, mean_ci, merge_levels, q_function)


# ------------------------------------------------------------------ analog rule
def gray_fields(data, side, per_alphabet, seed):
    """Score fields f = 1 - x of bilinear glyphs, with a fixed tiny jitter that makes them injective."""
    folder = data / "images_evaluation"
    rng = np.random.default_rng(seed)
    fields = []
    for alphabet in sorted(p for p in folder.iterdir() if p.is_dir() and not p.name.startswith(".")):
        paths = sorted(alphabet.glob("*/*.png"))
        for i in sorted(rng.choice(len(paths), per_alphabet, replace=False)):
            with Image.open(paths[i]) as src:
                ink = 1 - np.asarray(src.convert("L").resize((side, side), Image.Resampling.BILINEAR), float) / 255
            fields.append(0.01 + 0.98 * (1 - ink) + 1e-4 * rng.random((side, side)))  # Inside (0, 1), no ties.
    return fields


def relation_risk(levels, noisy_field):
    """Exact E |F(m_f) - F(m_fhat)| over ordered pairs for a uniform threshold on [0, 1]."""
    return float(np.abs(np.clip(levels, 0, 1) - np.clip(merge_levels(noisy_field), 0, 1)).mean())


def analog(data, sigmas, draws, seed):
    fields = gray_fields(data, 16, 2, seed)
    rng = np.random.default_rng(seed + 1)
    rows = {s: dict(uniform=[], allocated=[], linear_uniform=[], inside=[]) for s in sigmas}
    predicted = []
    for f in fields:
        load = critical_loads(f)
        power = load ** (2 / 3)
        power /= power.mean()
        levels = merge_levels(f)
        gap = np.diff(np.sort(f.ravel())).min()
        predicted.append(float(np.mean((load.size * load) ** (2 / 3)) ** 1.5))
        for s in sigmas:
            uni, alloc = [], []
            for _ in range(draws):
                noise = rng.standard_normal(f.shape)
                uni.append(relation_risk(levels, f + s * noise))
                alloc.append(relation_risk(levels, f + s * noise / np.sqrt(power)))
            rows[s]["uniform"].append(np.mean(uni))
            rows[s]["allocated"].append(np.mean(alloc))
            rows[s]["linear_uniform"].append(s * np.sqrt(2 / np.pi))   # sum_v c_v E|delta_v|
            rows[s]["inside"].append(3 * s < gap / 2)
    out = []
    for s in sigmas:
        r = rows[s]
        ratio = np.array(r["allocated"]) / np.array(r["uniform"])
        out.append(dict(sigma=s, snr_db=float(-20 * np.log10(s)), uniform_risk=float(np.mean(r["uniform"])),
                        allocated_risk=float(np.mean(r["allocated"])), ratio=mean_ci(ratio),
                        measured_over_linearized=float(np.mean(r["uniform"]) / np.mean(r["linear_uniform"])),
                        fields_within_half_gap=int(np.sum(r["inside"]))))
    return dict(fields=len(fields), side=16, draws=draws, predicted_ratio=float(np.mean(predicted)), rows=out)


# ------------------------------------------------------------------ shared power map
def training_glyphs(data, side, count, seed):
    folder = data / "images_background"
    paths = sorted(folder.glob("*/*/*.png"))
    rng = np.random.default_rng(seed)
    masks = []
    for i in sorted(rng.choice(len(paths), count, replace=False)):
        with Image.open(paths[i]) as src:
            masks.append(np.asarray(src.convert("L").resize((side, side), Image.Resampling.NEAREST)) < 128)
    return masks


def shared(data, snrs_db, draws, seed, train_count, test_per_alphabet):
    mean_influence = np.mean([influence_map(m) for m in training_glyphs(data, 32, train_count, seed)], axis=0)
    glyphs = load_glyphs(data, 32, test_per_alphabet, 20261115)
    maps = [influence_map(m) for _, m in glyphs]
    rng = np.random.default_rng(seed + 2)
    rows = []
    for snr_db in snrs_db:
        snr = 10 ** (snr_db / 10)
        common = bpsk_allocation(np.maximum(mean_influence, 1e-12), snr)
        risk = dict(uniform=[], per_image=[], shared=[], shuffled=[])
        errors = dict(uniform=[], per_image=[], shared=[], shuffled=[])
        for (_, mask), inf in zip(glyphs, maps):
            base = components(mask)
            # Same power values as per_image, assigned to randomly permuted locations.
            scrambled = rng.permutation(bpsk_allocation(inf, snr).ravel()).reshape(mask.shape)
            powers = dict(uniform=np.full(mask.shape, snr), per_image=bpsk_allocation(inf, snr), shared=common,
                          shuffled=scrambled)
            acc = {k: [] for k in powers}
            err = {k: [] for k in powers}
            for _ in range(draws):
                noise = rng.standard_normal(mask.shape)
                for name, p in powers.items():
                    flips = noise > np.sqrt(2 * p)
                    acc[name].append(disagreement(base, components(mask ^ flips)))
                    err[name].append(flips.mean())
            for name in powers:
                risk[name].append(np.mean(acc[name]))
                errors[name].append(np.mean(err[name]))
        u = np.mean(risk["uniform"])
        rows.append(dict(ebn0_db=snr_db, uniform_risk=float(u),
                         per_image_ratio=float(np.mean(risk["per_image"]) / u),
                         shared_ratio=float(np.mean(risk["shared"]) / u),
                         shuffled_ratio=float(np.mean(risk["shuffled"]) / u),
                         pixel_error={k: float(np.mean(v)) for k, v in errors.items()}))
    return dict(training_glyphs=train_count, test_glyphs=len(glyphs), draws=draws, rows=rows)


# ------------------------------------------------------------------ converse on orbits
def h2(p):
    p = np.clip(np.asarray(p, float), 1e-300, 1 - 1e-16)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))


def joint_entropy(answers):
    _, counts = np.unique(answers, axis=0, return_counts=True)
    p = counts / counts.sum()
    return float(-(p * np.log2(p)).sum())


def dictionary(labels, rng, d=3, candidates=4000):
    """Greedy choice of d location pairs whose answers best identify the state."""
    n = labels[0].size
    pairs = rng.integers(n, size=(candidates, 2))
    flat = np.stack([l.ravel() for l in labels])                       # [states, n]
    a, b = flat[:, pairs[:, 0]], flat[:, pairs[:, 1]]
    table = ((a > 0) & (a == b)).astype(np.int8)                       # [states, candidates]
    chosen = []
    for _ in range(d):
        best = max(range(candidates), key=lambda j: joint_entropy(table[:, chosen + [j]]))
        chosen.append(best)
    return table[:, chosen]


def bpsk_capacity(snr, nodes=4001):
    """Mutual information in bits of BPSK on a complex AWGN use with symbol SNR `snr`."""
    a = np.sqrt(2 * snr)
    y = np.linspace(-12, 12, nodes)
    f = np.exp(-(y - a) ** 2 / 2) / np.sqrt(2 * np.pi)
    return float(1 - np.sum(f * np.log2(1 + np.exp(-2 * a * y))) * (y[1] - y[0]))


def converse(data, snrs_db, seed):
    glyphs = load_glyphs(data, 32, 16, 20261115)
    rng = np.random.default_rng(seed + 3)
    tables = []
    for _, mask in glyphs:
        members = orbit(mask)
        if len(members) == 8:
            tables.append(dictionary([components(m) for m in members], rng))
    hq = np.array([joint_entropy(t) for t in tables])
    rows = []
    for snr_db in snrs_db:
        snr = 10 ** (snr_db / 10)
        bit_error = float(q_function(np.sqrt(2 * snr)))
        # Three BPSK symbols carry the state index; a received index is any of the 8 states.
        index = np.arange(8)
        hamming = np.array([[bin(i ^ j).count("1") for j in index] for i in index])
        transition = bit_error**hamming * (1 - bit_error) ** (3 - hamming)
        deltas = []
        for t in tables:
            errors = np.array([np.mean(np.sum(transition * (t[:, None, j] != t[None, :, j]), axis=1))
                               for j in range(t.shape[1])])
            deltas.append(joint_entropy(t) - h2(errors).sum())
        need = np.maximum(deltas, 0)
        rows.append(dict(snr_db=snr_db, bit_error=bit_error, mean_bound_bits=float(need.mean()),
                         information_sent_bits=3 * bpsk_capacity(snr), capacity_bits=3 * float(np.log2(1 + snr)),
                         min_uses=float((need / np.log2(1 + snr)).mean()), uses=3,
                         violations=int(np.sum(need > 3 * bpsk_capacity(snr) + 1e-9))))
    return dict(orbits=len(tables), queries=3, mean_answer_entropy_bits=float(hq.mean()),
                share_with_three_bits=float(np.mean(hq > 3 - 1e-9)), rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--only", choices=("analog", "shared", "converse"))
    args = parser.parse_args()
    record = {}
    if args.only in (None, "converse"):
        record["converse"] = converse(args.data, (0, 3, 6, 10), args.seed)
        print(json.dumps(record["converse"], indent=1), flush=True)
    if args.only in (None, "analog"):
        record["analog"] = analog(args.data, (1e-6, 1e-3, 0.01, 0.03, 0.1, 0.3), 20, args.seed)
        print(json.dumps(record["analog"], indent=1), flush=True)
    if args.only in (None, "shared"):
        record["shared"] = shared(args.data, (4, 6, 8), 100, args.seed, 200, 16)
        print(json.dumps(record["shared"], indent=1), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=1) + "\n")
