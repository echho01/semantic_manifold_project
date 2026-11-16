"""Relation rate-distortion on a symmetry orbit of the readout grid.

A template and its images under grid automorphisms share every unlabelled
persistence diagram. This script evaluates, for binary masks, the exact
pair-disagreement matrix on the dihedral orbit and the closed-form orbit RDF
R_Q(delta) = log2 K - H_2(P_beta). It also tabulates the channel-use converse
for complex AWGN and ergodic coherent Rayleigh fading (formula evaluation).
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import exp1

from reflection_queries import discordance, foreground_components

# The eight automorphisms of a square four-neighbor grid, as array maps.
DIHEDRAL = {
    "id": lambda a: a,
    "rot90": lambda a: np.rot90(a, 1),
    "rot180": lambda a: np.rot90(a, 2),
    "rot270": lambda a: np.rot90(a, 3),
    "flip_lr": np.fliplr,
    "flip_ud": np.flipud,
    "transpose": lambda a: a.T,
    "anti_transpose": lambda a: np.rot90(a, 2).T,
}


def orbit(mask, names=tuple(DIHEDRAL)):
    """Distinct images of a mask under the named automorphisms, in first-seen order."""
    seen, members = set(), []
    for name in names:
        image = np.ascontiguousarray(DIHEDRAL[name](mask))
        key = image.tobytes()
        if key not in seen:
            seen.add(key)
            members.append(image)
    return members


def distortion_matrix(members):
    """Exact uniform-pair, uniform-threshold disagreement between binary masks."""
    labels = [foreground_components(m) for m in members]
    eligible = np.ones(members[0].shape, bool)
    k = len(members)
    rho = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            rho[i, j] = rho[j, i] = discordance(labels[i], labels[j], eligible)
    return rho


def gibbs(row, beta):
    w = np.exp2(-beta * (row - row.min()))
    return w / w.sum()


def entropy_bits(p):
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def orbit_rate(row, delta):
    """log2 K - H_2(P_beta) with the Gibbs mean equal to delta (bits)."""
    row = np.asarray(row, float)
    k = len(row)
    if delta >= row.mean():
        return 0.0
    zeros = np.isclose(row, 0.0)
    if delta <= 0:
        return float(np.log2(k) - np.log2(zeros.sum()))
    lo, hi = 0.0, 1.0
    while gibbs(row, hi) @ row > delta:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if gibbs(row, mid) @ row > delta:
            lo = mid
        else:
            hi = mid
    return float(np.log2(k) - entropy_bits(gibbs(row, hi)))


def blahut_point(rho, beta, iters=20000, tol=1e-13):
    """(distortion, rate in bits) of the Blahut optimum at slope beta, uniform source."""
    k = len(rho)
    p = np.full(k, 1 / k)
    r = np.full(k, 1 / k)
    kernel = np.exp2(-beta * rho)
    for _ in range(iters):
        w = r * kernel
        w /= w.sum(1, keepdims=True)
        new = p @ w
        if np.abs(new - r).max() < tol:
            r = new
            break
        r = new
    w = r * kernel
    w /= w.sum(1, keepdims=True)
    joint = p[:, None] * w
    marg = joint.sum(0)
    mask = joint > 0
    rate = float((joint[mask] * np.log2((w / marg)[mask])).sum())
    return float((joint * rho).sum()), rate


def channel_table(snrs_db):
    rows = []
    for snr in snrs_db:
        lin = 10 ** (snr / 10)
        awgn = np.log2(1 + lin)
        rayleigh = float(np.exp(1 / lin) * exp1(1 / lin) / np.log(2))
        rows.append(dict(snr_db=snr, awgn_bits_per_use=float(awgn),
                         rayleigh_ergodic_bits_per_use=rayleigh,
                         awgn_uses_per_bit=float(1 / awgn),
                         rayleigh_uses_per_bit=1 / rayleigh))
    return rows


def synthetic(seed, side, count, density):
    """Random binary templates; exact formula against Blahut on each orbit."""
    rng = np.random.default_rng(seed)
    records, worst = [], 0.0
    while len(records) < count:
        mask = rng.random((side, side)) < density
        members = orbit(mask)
        rho = distortion_matrix(members)
        k = len(members)
        sorted_rows = np.sort(rho, axis=1)
        balanced = bool(np.allclose(sorted_rows, sorted_rows[0]))
        gaps = []
        for beta in (2.0, 8.0, 32.0, 128.0):
            d, r = blahut_point(rho, beta)
            gaps.append(abs(r - orbit_rate(rho[0], d)))
        worst = max(worst, max(gaps))
        records.append(dict(orbit_size=k, balanced_rows=balanced,
                            delta_max=float(rho[0].mean()),
                            max_formula_minus_blahut=float(max(gaps))))
    return dict(seed=seed, side=side, count=count, density=density,
                all_balanced=all(r["balanced_rows"] for r in records),
                worst_formula_minus_blahut_bits=worst, templates=records)


GRID = np.linspace(0.0, 1.0, 41)


def plot_curves(curve, output):
    import matplotlib.pyplot as plt

    grid = np.array(curve["fraction"])
    fig, ax = plt.subplots(figsize=(3.4, 2.1), layout="constrained")
    ax.fill_between(grid, curve["p10"], curve["p90"], color="#123f63", alpha=0.2, linewidth=0)
    ax.plot(grid, curve["mean"], color="#123f63", label="Dihedral orbit, $K=8$")
    ax.plot(grid, curve["pair"], color="#a32437", linestyle="--", label="Reflection pair, $K=2$")
    ax.axhline(0.0, color="0.35", linewidth=0.8, linestyle=":", label="Any diagram distortion")
    ax.set(xlabel=r"Relation distortion $\delta/\delta_{\max}$", ylabel="Rate (bit/glyph)",
           xlim=(0, 1), ylim=(-0.1, 3.1))
    ax.legend(frameon=False, fontsize=7)
    ax.tick_params(labelsize=7)
    for label in (ax.xaxis.label, ax.yaxis.label):
        label.set_fontsize(8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def omniglot(data, side, per_alphabet, seed, reference):
    from PIL import Image

    folder = data / "images_evaluation"
    alphabets = sorted(p for p in folder.iterdir() if p.is_dir() and not p.name.startswith("."))
    rng = np.random.default_rng(seed)
    records, curves = [], []
    for alphabet in alphabets:
        paths = sorted(alphabet.glob("*/*.png"))
        for i in sorted(rng.choice(len(paths), per_alphabet, replace=False)):
            path = paths[i]
            with Image.open(path) as src:
                mask = np.asarray(src.convert("L").resize((side, side),
                                                          Image.Resampling.NEAREST)) < 128
            members = orbit(mask)
            rho = distortion_matrix(members)
            row = rho[0]
            if len(members) == 8:
                curves.append([orbit_rate(row, f * row.mean()) for f in GRID])
            flip = orbit(mask, ("id", "flip_lr"))
            alpha = float(distortion_matrix(flip)[0, 1]) if len(flip) == 2 else 0.0
            dmax = float(row.mean())
            records.append(dict(
                id=path.relative_to(folder).as_posix(), orbit_size=len(members),
                reflection_alpha=alpha, delta_max=dmax,
                min_offdiag=float(row[1:].min()) if len(row) > 1 else None,
                max_offdiag=float(row[1:].max()) if len(row) > 1 else None,
                rate_bits={str(f): orbit_rate(row, f * dmax) for f in (0.0, 0.1, 0.25, 0.5)},
                reflection_rate_bits={str(f): orbit_rate(np.array([0.0, alpha]), f * alpha / 2)
                                      for f in (0.0, 0.1, 0.25, 0.5)} if alpha > 0 else None))
    ids = hashlib.sha256("\n".join(r["id"] for r in records).encode()).hexdigest()
    alphas = np.array([r["reflection_alpha"] for r in records])
    summary = dict(seed=seed, side=side, count=len(records), split="official evaluation",
                   id_sha256=ids, reflection_alpha_mean=float(alphas.mean()),
                   reflection_alpha_median=float(np.median(alphas)))
    if reference:
        ref = json.loads(reference.read_text("utf-8"))
        summary["reference_id_sha256_match"] = ref["id_sha256"] == ids
        by_id = {r["id"]: r["uniform"] for r in ref["images"]}
        summary["reference_max_abs_alpha_difference"] = float(max(
            abs(by_id[r["id"]] - r["reflection_alpha"]) for r in records))
    sizes = np.array([r["orbit_size"] for r in records])
    summary["orbit_size_counts"] = {str(int(s)): int((sizes == s).sum()) for s in np.unique(sizes)}
    full = [r for r in records if r["orbit_size"] == 8]
    dmax = np.array([r["delta_max"] for r in full])
    summary["full_orbit"] = dict(
        count=len(full), delta_max_mean=float(dmax.mean()), delta_max_median=float(np.median(dmax)),
        min_offdiag_mean=float(np.mean([r["min_offdiag"] for r in full])),
        max_offdiag_mean=float(np.mean([r["max_offdiag"] for r in full])),
        rate_bits_mean={f: float(np.mean([r["rate_bits"][f] for r in full]))
                        for f in ("0.0", "0.1", "0.25", "0.5")},
        rate_bits_median={f: float(np.median([r["rate_bits"][f] for r in full]))
                          for f in ("0.0", "0.1", "0.25", "0.5")})
    pair = [r for r in records if r["reflection_rate_bits"]]
    summary["reflection_pair"] = dict(
        count=len(pair),
        rate_bits_mean={f: float(np.mean([r["reflection_rate_bits"][f] for r in pair]))
                        for f in ("0.0", "0.1", "0.25", "0.5")})
    curves = np.array(curves)
    summary["curve"] = dict(
        fraction=GRID.tolist(), mean=curves.mean(0).tolist(),
        p10=np.quantile(curves, 0.1, axis=0).tolist(), p90=np.quantile(curves, 0.9, axis=0).tolist(),
        pair=[orbit_rate(np.array([0.0, 1.0]), f / 2) for f in GRID])
    return dict(**summary, images=records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--data", type=Path, help="Omniglot root; omit for synthetic checks only")
    ap.add_argument("--reference", type=Path, help="reflection_queries.json for a pipeline check")
    ap.add_argument("--figure", type=Path, help="write the orbit rate curves (needs --data)")
    ap.add_argument("--side", type=int, default=32)
    ap.add_argument("--per-alphabet", type=int, default=16)
    ap.add_argument("--seed", type=int, default=20261115)
    ap.add_argument("--synthetic-side", type=int, default=5)
    ap.add_argument("--synthetic-count", type=int, default=200)
    a = ap.parse_args()
    result = dict(
        synthetic=synthetic(a.seed, a.synthetic_side, a.synthetic_count, 0.45),
        channel=channel_table([0.0, 5.0, 10.0, 15.0, 20.0]))
    if a.data:
        result["omniglot"] = omniglot(a.data, a.side, a.per_alphabet, a.seed, a.reference)
        if a.figure:
            plot_curves(result["omniglot"]["curve"], a.figure)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + "\n", "utf-8")
    brief = {k: ({kk: vv for kk, vv in v.items() if kk not in ("images", "templates", "curve")}
                 if isinstance(v, dict) else v) for k, v in result.items()}
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
