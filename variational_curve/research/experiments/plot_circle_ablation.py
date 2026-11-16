"""Plot achieved rate and the excess over the finite BA curve at matched distortion."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    research = Path(__file__).resolve().parents[1]
    with (research / "code/results/early_surface/G5_geometry_density.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    costs = np.linspace(0.1, 1.0, 241)
    def rate(mode):
        points = sorted((float(r[mode + "_D"]), float(r[mode + "_R"])) for r in rows)
        d, bits = np.asarray(points).T
        if costs.min() < d.min() or costs.max() > d.max():
            raise ValueError("matched distortion outside the recorded range")
        return np.interp(costs, d, bits)
    ref = rate("BA")
    fig, (left, right) = plt.subplots(1, 2, figsize=(6.8, 2.35), layout="constrained")
    styles = [("fixed", "#30495e", "-", "Fixed"), ("geometry", "#b75d16", "--", "Geometry"),
              ("density", "#267260", "-", "Density"), ("joint", "#782859", ":", "Joint")]
    for mode, color, linestyle, label in styles:
        values = rate(mode)
        left.plot(costs, values, color=color, ls=linestyle, label=label, lw=1.6)
        right.plot(costs, values - ref, color=color, ls=linestyle, lw=1.6)
    left.set(ylabel="Information (bit/state)", title="(a) Achievable information")
    right.set(ylabel="Rate minus BA (bit/state)", title="(b) Excess at matched distortion")
    for axis in (left, right):
        axis.axvline(0.30, color=".65", ls="--", lw=.8, zorder=0)
        axis.set(xlabel="Squared arc distortion", xlim=(.1, 1))
        axis.grid(alpha=.2)
        axis.tick_params(labelsize=7)
        axis.title.set_fontsize(9)
        axis.xaxis.label.set_fontsize(8)
        axis.yaxis.label.set_fontsize(8)
    left.legend(fontsize=7, ncol=2, frameon=False)
    dest = research.parent / "figures/circle_ablation.png"
    dest.parent.mkdir(exist_ok=True)
    fig.savefig(dest, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
