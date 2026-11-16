"""Exact connectivity disagreement under horizontal reflection of Omniglot masks."""

import argparse
import hashlib
import json
from collections import deque
from pathlib import Path

import numpy as np
from PIL import Image


def foreground_components(mask):
    height, width = mask.shape
    labels = np.zeros(mask.shape, dtype=np.int32)
    component = 0
    for row, col in zip(*np.nonzero(mask)):
        if labels[row, col]:
            continue
        component += 1
        queue = deque([(row, col)])
        labels[row, col] = component
        while queue:
            y, x = queue.popleft()
            for a, b in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
                if 0 <= a < height and 0 <= b < width and mask[a, b] and not labels[a, b]:
                    labels[a, b] = component
                    queue.append((a, b))
    return labels


def discordance(first, second, eligible):
    """Count ordered pairs with replacement. Zero labels are always inactive."""
    a, b = first[eligible], second[eligible]
    if not a.size:
        return None
    def square_sum(labels):
        _, counts = np.unique(labels, return_counts=True)
        return int(counts @ counts)
    active = (a > 0) & (b > 0)
    stride = int(b.max()) + 1
    intersection = square_sum(a[active].astype(np.int64) * stride + b[active])
    mismatches = square_sum(a[a > 0]) + square_sum(b[b > 0]) - 2 * intersection
    return mismatches / a.size**2


def verify():
    labels = foreground_components(np.array([[1, 0, 1], [1, 0, 0]], dtype=bool))
    flipped = np.fliplr(labels)
    for eligible in (np.ones(labels.shape, bool), labels > 0):
        a, b = labels[eligible], flipped[eligible]
        def connections(x):
            return (x[:, None] == x) & (x[:, None] > 0)
        assert np.isclose(discordance(labels, flipped, eligible),
                          np.mean(connections(a) != connections(b)))
    assert discordance(labels, labels, labels > 0) == 0.0


def plot_diagnostic(folder, records, side, output):
    import matplotlib.pyplot as plt

    values = [r for r in records if r["foreground"] is not None]
    median = np.median([r["foreground"] for r in values])
    # Prefer witnesses whose two endpoints stay on ink under reflection.
    for row in sorted(values, key=lambda r: abs(r["foreground"] - median)):
        with Image.open(folder / row["id"]) as image:
            mask = np.asarray(image.convert("L").resize((side, side), Image.Resampling.NEAREST)) < 128
        left = foreground_components(mask)
        right = np.fliplr(left)
        shared = np.argwhere((left > 0) & (right > 0))
        if not len(shared):
            continue
        a, b = left[tuple(shared.T)], right[tuple(shared.T)]
        candidates = np.argwhere((a[:, None] == a) & (b[:, None] != b))
        if candidates.size:
            distances = np.sum((shared[candidates[:, 0]] - shared[candidates[:, 1]])**2, axis=1)
            endpoints = shared[candidates[distances.argmax()]]
            break
    else:
        raise ValueError("no foreground pair witnesses a component change")
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.0),
                             gridspec_kw={"width_ratios": [1, 1, 1.65]}, layout="constrained")
    for ax, glyph, title in zip(axes[:2], (mask, np.fliplr(mask)),
                                ("Original: connected", "Reflected: disconnected")):
        ax.imshow(glyph, cmap="Greys", interpolation="nearest", vmin=0, vmax=1)
        ax.scatter(endpoints[:, 1], endpoints[:, 0], s=24, facecolors="white",
                   edgecolors="#a32437", linewidths=1.4)
        ax.set(title=title, xticks=[], yticks=[])
        ax.title.set_fontsize(8.5)
    ax = axes[2]
    for key, label, color in (("uniform", "All pixels", "#123f63"),
                              ("foreground", "Source foreground", "#a32437")):
        vals = sorted(r[key] for r in values)
        ax.step(vals, np.arange(1, len(vals) + 1) / len(vals), where="post", label=label, color=color)
    ax.set(xlabel="Disagreement probability", ylabel="Image fraction", xlim=(0, 1), ylim=(0, 1.02))
    ax.legend(frameon=False, fontsize=7, loc="center right")
    ax.tick_params(labelsize=7)
    for label in (ax.xaxis.label, ax.yaxis.label):
        label.set_fontsize(8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--per-alphabet", type=int, default=16)
    parser.add_argument("--side", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261115)
    args = parser.parse_args()
    if min(args.side, args.per_alphabet) < 1:
        parser.error("side and per-alphabet must be positive")
    verify()
    folder = args.data / "images_evaluation"
    alphabets = sorted(p for p in folder.iterdir() if p.is_dir() and not p.name.startswith("."))
    if not alphabets:
        parser.error("missing Omniglot evaluation alphabets")
    rng = np.random.default_rng(args.seed)
    records = []
    for alphabet in alphabets:
        paths = sorted(alphabet.glob("*/*.png"))
        for i in sorted(rng.choice(len(paths), args.per_alphabet, replace=False)):
            path = paths[i]
            with Image.open(path) as src:
                mask = np.asarray(src.convert("L").resize((args.side, args.side),
                                                          Image.Resampling.NEAREST)) < 128
            labels = foreground_components(mask)
            reflected = np.fliplr(labels)
            uniform = discordance(labels, reflected, np.ones_like(mask))
            foreground = discordance(labels, reflected, mask)
            records.append(dict(id=path.relative_to(folder).as_posix(),
                                uniform=uniform, foreground=foreground,
                                active_fraction=float(mask.mean())))
    result = dict(seed=args.seed, side=args.side, count=len(records), alphabets=len(alphabets),
                  per_alphabet=args.per_alphabet, split="official evaluation",
                  id_sha256=hashlib.sha256("\n".join(r["id"] for r in records).encode()).hexdigest(),
                  nonempty=sum(r["foreground"] is not None for r in records))
    for name in ("uniform", "foreground"):
        risks = np.array([r[name] for r in records if r[name] is not None])
        result[name] = dict(mean=float(risks.mean()), median=float(np.median(risks)),
                            positive=int(np.count_nonzero(risks)), count=len(risks))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(**result, images=records), indent=2) + "\n", "utf-8")
    if args.figure:
        plot_diagnostic(folder, records, args.side, args.figure)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
