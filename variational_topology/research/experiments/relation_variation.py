"""First variation of the relation functional: critical loads, influences, allocations.

For a score field f on the four-neighbor grid, the relation functional is the
pair-disagreement d_nu(f, f_hat) under uniform ordered location pairs (with
replacement) and a uniform threshold. This script evaluates

  * critical loads of generic (injective) fields by one union-find sweep;
  * exact single-pixel influences Inf_x(v) = d_nu(x, x + e_v) of binary glyphs;
  * the first-order law E d_nu(x, x + Z) ~ sum_v eps_v Inf_x(v) by Monte Carlo;
  * an influence-weighted BPSK power allocation against uniform power, with
    common noise for the two schemes.

The Omniglot protocol is that of reflection_queries.py and orbit_rdf.py.
The link experiment is uncoded and controlled. It does not rank neural codecs.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.special import erfc, lambertw

from reflection_queries import discordance, foreground_components


# ---------------------------------------------------------------- generic fields
def grid_neighbors(height, width):
    idx = np.arange(height * width).reshape(height, width)
    nbrs = [[] for _ in range(height * width)]
    for a, b in ((idx[:-1], idx[1:]), (idx[:, :-1], idx[:, 1:])):
        for u, v in zip(a.ravel(), b.ravel()):
            nbrs[u].append(int(v))
            nbrs[v].append(int(u))
    return nbrs


def critical_loads(field):
    """Loads c_f(v) of an injective field under uniform ordered pairs with replacement.

    Vertices enter in increasing order of f. A vertex that joins components of
    sizes n_1, n_2, ... first connects N^2 - sum n_j^2 ordered pairs, N = 1 + sum n_j.
    """
    f = np.asarray(field, float)
    flat = f.ravel()
    if len(np.unique(flat)) != flat.size:
        raise ValueError("critical loads need an injective field")
    nbrs = grid_neighbors(*f.shape)
    parent = np.full(flat.size, -1)
    size = np.zeros(flat.size, np.int64)
    load = np.zeros(flat.size)

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for v in np.argsort(flat):
        roots = {find(u) for u in nbrs[v] if parent[u] >= 0}
        sizes = [int(size[r]) for r in roots]
        total = 1 + sum(sizes)
        load[v] = total**2 - sum(s * s for s in sizes)
        parent[v] = v
        size[v] = total
        for r in roots:
            parent[r] = v
    return (load / flat.size**2).reshape(f.shape)


def merge_levels(field):
    """All-pairs minimax levels m_f(a, b) by closure; for brute-force checks only."""
    f = np.asarray(field, float)
    flat = f.ravel()
    n = flat.size
    m = np.full((n, n), np.inf)
    m[np.arange(n), np.arange(n)] = flat
    for u, vs in enumerate(grid_neighbors(*f.shape)):
        for v in vs:
            m[u, v] = max(flat[u], flat[v])
    for k in range(n):
        m = np.minimum(m, np.maximum(m[:, k:k + 1], m[k:k + 1, :]))
    return m


def analog_gain(load):
    """Risk ratio of the p ~ c^(2/3) allocation to uniform power, Gaussian perturbation."""
    c = np.asarray(load, float).ravel()
    return float(np.mean((c.size * c) ** (2 / 3)) ** 1.5)


# ---------------------------------------------------------------- binary images
FOUR = ndimage.generate_binary_structure(2, 1)


def components(mask):
    return ndimage.label(mask, structure=FOUR)[0]


def disagreement(first, second):
    """Exact d_nu between two label images, all pixels eligible (zero = inactive)."""
    a, b = first.ravel().astype(np.int64), second.ravel().astype(np.int64)
    ca, cb = np.bincount(a)[1:], np.bincount(b)[1:]
    both = (a > 0) & (b > 0)
    joint = np.bincount(a[both] * (int(b.max()) + 1) + b[both])
    return float(ca @ ca + cb @ cb - 2 * (joint @ joint)) / a.size**2


def influence_map(mask):
    """Inf_x(v) = d_nu(x, x + e_v) for every pixel v."""
    base = components(mask)
    out = np.zeros(mask.shape)
    work = mask.copy()
    for v in np.ndindex(mask.shape):
        work[v] = not work[v]
        out[v] = disagreement(base, components(work))
        work[v] = not work[v]
    return out


def bpsk_allocation(weights, mean_snr):
    """Minimize sum_v w_v Q(sqrt(2 p_v)) subject to mean(p) = mean_snr (N0 = 1).

    Stationarity gives p_v = W(w_v^2 / (2 pi mu^2)) / 2 with the Lambert function W.
    """
    w = np.asarray(weights, float).ravel()
    if np.any(w <= 0):
        raise ValueError("every pixel has positive influence")

    def powers(log_mu):
        return 0.5 * np.real(lambertw(np.exp(2 * (np.log(w) - log_mu)) / (2 * np.pi)))

    lo, hi = np.log(w.min()) - 60.0, np.log(w.max()) + 60.0     # mean power decreases with mu
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if powers(mid).mean() > mean_snr:
            lo = mid
        else:
            hi = mid
    return powers(0.5 * (lo + hi)).reshape(np.shape(weights))


def q_function(x):
    return 0.5 * erfc(np.asarray(x) / np.sqrt(2))


# ---------------------------------------------------------------- Omniglot study
def load_glyphs(data, side, per_alphabet, seed):
    from PIL import Image

    folder = data / "images_evaluation"
    alphabets = sorted(p for p in folder.iterdir() if p.is_dir() and not p.name.startswith("."))
    rng = np.random.default_rng(seed)
    glyphs = []
    for alphabet in alphabets:
        paths = sorted(alphabet.glob("*/*.png"))
        for i in sorted(rng.choice(len(paths), per_alphabet, replace=False)):
            with Image.open(paths[i]) as src:
                mask = np.asarray(src.convert("L").resize((side, side), Image.Resampling.NEAREST)) < 128
            glyphs.append((paths[i].relative_to(folder).as_posix(), mask))
    return glyphs


def mean_ci(values):
    v = np.asarray(values, float)
    half = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
    return dict(mean=float(v.mean()), ci95=[float(v.mean() - half), float(v.mean() + half)])


def omniglot(data, side, per_alphabet, seed, reference, flip_draws, link_draws, snrs_db, maps_out, figure):
    glyphs = load_glyphs(data, side, per_alphabet, seed)
    ids = hashlib.sha256("\n".join(g[0] for g in glyphs).encode()).hexdigest()
    n = side * side

    # Pipeline check against the stored reflection experiment.
    all_px, fg = [], []
    for _, mask in glyphs:
        labels = foreground_components(mask)
        if not np.isclose(disagreement(labels, np.fliplr(labels)),
                          discordance(labels, np.fliplr(labels), np.ones(mask.shape, bool))):
            raise ArithmeticError("fast pair count disagrees with the reference implementation")
        all_px.append(disagreement(labels, np.fliplr(labels)))
        fg.append(discordance(labels, np.fliplr(labels), labels > 0))
    check = dict(id_sha256=ids, reflection_all_pixels_mean=float(np.mean(all_px)),
                 reflection_foreground_mean=float(np.mean(fg)))
    if reference:
        ref = json.loads(reference.read_text("utf-8"))
        by_id = {r["id"]: r for r in ref["images"]}
        check["reference_id_sha256_match"] = ref["id_sha256"] == ids
        check["reference_max_abs_difference"] = float(max(
            max(abs(by_id[g[0]]["uniform"] - a), abs(by_id[g[0]]["foreground"] - f))
            for g, a, f in zip(glyphs, all_px, fg)))

    # Exact influence maps.
    maps = np.stack([influence_map(mask) for _, mask in glyphs])
    masks = np.stack([mask for _, mask in glyphs])
    records = []
    for (name, mask), inf in zip(glyphs, maps):
        total = float(inf.sum())
        ranked = np.sort(inf.ravel())[::-1]
        records.append(dict(
            id=name, foreground_pixels=int(mask.sum()), total_influence=total,
            foreground_influence=float(inf[mask].sum()), background_influence=float(inf[~mask].sum()),
            max_influence=float(ranked[0]), min_influence=float(ranked[-1]),
            top1_share=float(ranked[:round(0.01 * n)].sum() / total),
            top5_share=float(ranked[:round(0.05 * n)].sum() / total)))
    influence = dict(
        pixels=n, top1_pixels=round(0.01 * n), top5_pixels=round(0.05 * n),
        total=mean_ci([r["total_influence"] for r in records]),
        total_median=float(np.median([r["total_influence"] for r in records])),
        foreground_share_mean=float(np.mean([r["foreground_influence"] / r["total_influence"] for r in records])),
        foreground_pixel_share_mean=float(np.mean([r["foreground_pixels"] / n for r in records])),
        top1_share_mean=float(np.mean([r["top1_share"] for r in records])),
        top5_share_mean=float(np.mean([r["top5_share"] for r in records])),
        max_influence_mean=float(np.mean([r["max_influence"] for r in records])),
        uniform_share_top1=round(0.01 * n) / n, uniform_share_top5=round(0.05 * n) / n)
    if maps_out:
        maps_out.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(maps_out, ids=np.array([g[0] for g in glyphs]), masks=masks, influence=maps)

    rng = np.random.default_rng(seed + 1)
    bases = [components(mask) for mask in masks]

    # First-order law under independent flips of equal probability.
    first_order = []
    for eps, draws in flip_draws:
        simulated = []
        for mask, base in zip(masks, bases):
            flips = rng.random((draws,) + mask.shape) < eps
            risks = [disagreement(base, components(mask ^ z)) if z.any() else 0.0 for z in flips]
            simulated.append(np.mean(risks))
        predicted = eps * np.array([r["total_influence"] for r in records])
        ratio = float(np.mean(simulated) / predicted.mean())
        first_order.append(dict(epsilon=eps, draws_per_glyph=draws, expected_flips=eps * n,
                                worst_case_remainder=(eps * n) ** 2,
                                predicted_mean=float(predicted.mean()), simulated=mean_ci(simulated),
                                paired_simulated_minus_predicted=mean_ci(np.array(simulated) - predicted),
                                simulated_over_predicted=ratio))

    # Influence-weighted BPSK power against uniform power, common noise.
    link = []
    for snr_db in snrs_db:
        snr = 10 ** (snr_db / 10)
        uni_risk, var_risk, uni_ber, var_ber, pred_uni, pred_var, spread = [], [], [], [], [], [], []
        for mask, base, inf in zip(masks, bases, maps):
            p = bpsk_allocation(inf, snr)
            noise = rng.standard_normal((link_draws,) + mask.shape)
            eu = noise > np.sqrt(2 * snr)
            ev = noise > np.sqrt(2 * p)
            uni_risk.append(np.mean([disagreement(base, components(mask ^ z)) if z.any() else 0.0 for z in eu]))
            var_risk.append(np.mean([disagreement(base, components(mask ^ z)) if z.any() else 0.0 for z in ev]))
            uni_ber.append(eu.mean())
            var_ber.append(ev.mean())
            pred_uni.append(float(inf.sum() * q_function(np.sqrt(2 * snr))))
            pred_var.append(float((inf * q_function(np.sqrt(2 * p))).sum()))
            spread.append(float(10 * np.log10(p.max() / p.min())))
        diff = np.array(uni_risk) - np.array(var_risk)
        link.append(dict(
            ebn0_db=snr_db, draws_per_glyph=link_draws,
            uniform_risk=mean_ci(uni_risk), allocated_risk=mean_ci(var_risk),
            paired_reduction=mean_ci(diff), risk_ratio=float(np.mean(var_risk) / np.mean(uni_risk)),
            first_order_uniform=float(np.mean(pred_uni)), first_order_allocated=float(np.mean(pred_var)),
            uniform_pixel_error=float(np.mean(uni_ber)), allocated_pixel_error=float(np.mean(var_ber)),
            theory_uniform_pixel_error=float(q_function(np.sqrt(2 * snr))),
            power_spread_db_mean=float(np.mean(spread))))

    shown = plot_influence(glyphs, maps, records, figure) if figure else None
    return dict(figure_glyph=shown, seed=seed, side=side, count=len(glyphs), split="official evaluation",
                pipeline_check=check, influence=influence, first_order=first_order, bpsk=link,
                images=records)


def plot_influence(glyphs, maps, records, output):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    totals = np.array([r["total_influence"] for r in records])
    pick = int(np.argsort(totals)[len(totals) // 2])                 # glyph of median total influence
    mask, inf = glyphs[pick][1], maps[pick]
    fig, axes = plt.subplots(1, 2, figsize=(3.4, 1.75), layout="constrained")
    axes[0].imshow(mask, cmap="gray_r", interpolation="nearest")
    axes[0].set_title("Glyph", fontsize=8)
    im = axes[1].imshow(inf, norm=LogNorm(vmin=inf.min(), vmax=inf.max()), cmap="magma",
                        interpolation="nearest")
    axes[1].set_title(r"Influence $\mathrm{Inf}_x(v)$", fontsize=8)
    bar = fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    bar.ax.tick_params(labelsize=6)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return glyphs[pick][0]


# ---------------------------------------------------------------- synthetic checks
def synthetic(seed, side, count):
    """Brute-force checks of the load formula and the local weighted-l1 law."""
    rng = np.random.default_rng(seed)
    worst_law = worst_sum = worst_count = 0.0
    gains = []
    for _ in range(count):
        f = rng.permutation(side * side).reshape(side, side) / (side * side) + 0.25 / (side * side)
        gap = 1.0 / (side * side)
        load = critical_loads(f)
        m = merge_levels(f)
        # First connecting vertex of each pair, from the brute-force levels.
        owner = np.searchsorted(np.sort(f.ravel()), m)
        brute = np.bincount(np.argsort(f.ravel())[owner].ravel(), minlength=f.size) / f.size**2
        worst_count = max(worst_count, float(np.abs(brute - load.ravel()).max()))
        worst_sum = max(worst_sum, abs(float(load.sum()) - 1.0))
        delta = rng.uniform(-0.49 * gap, 0.49 * gap, f.shape)
        direct = float(np.abs(merge_levels(f + delta) - m).mean())
        worst_law = max(worst_law, abs(direct - float((load * np.abs(delta)).sum())))
        gains.append(analog_gain(load))
    return dict(seed=seed, side=side, count=count, max_load_formula_error=worst_count,
                max_load_sum_error=worst_sum, max_local_law_error=worst_law,
                analog_gain_mean=float(np.mean(gains)), analog_gain_min=float(np.min(gains)),
                analog_gain_max=float(np.max(gains)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--data", type=Path, help="Omniglot root; omit for synthetic checks only")
    ap.add_argument("--reference", type=Path, help="reflection_queries.json for a pipeline check")
    ap.add_argument("--maps", type=Path, help="write per-glyph influence maps (.npz)")
    ap.add_argument("--figure", type=Path, help="write a glyph with its influence map")
    ap.add_argument("--side", type=int, default=32)
    ap.add_argument("--per-alphabet", type=int, default=16)
    ap.add_argument("--seed", type=int, default=20261115)
    ap.add_argument("--link-draws", type=int, default=100)
    a = ap.parse_args()
    result = dict(synthetic=synthetic(a.seed, 6, 100))
    if a.data:
        result["omniglot"] = omniglot(
            a.data, a.side, a.per_alphabet, a.seed, a.reference,
            flip_draws=((1e-4, 1000), (1e-3, 300), (1e-2, 100)), link_draws=a.link_draws,
            snrs_db=(4.0, 6.0, 8.0), maps_out=a.maps, figure=a.figure)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + "\n", "utf-8")
    brief = {k: ({kk: vv for kk, vv in v.items() if kk != "images"} if isinstance(v, dict) else v)
             for k, v in result.items()}
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
