"""Knowledge-conditioned rate-distortion on frozen CIFAR-100 representations.

The knowledge base is K = (theta, S).
  theta  parametric knowledge: the frozen backbone and the task head. It is the
         same at both terminals for every sample, so it carries no I(U; theta).
         It fixes the representation cloud, the cost and the anchors.
  S      situational knowledge: a discrete variable that changes per sample and
         is given to both terminals. Here it is the dataset superclass (20
         values) or the fine class (100 values). The dataset label is an oracle
         stand-in for context that a deployed pair would share. It is never sent.

For each theta the same source states, the same anchors and the same cost are
used with S absent, S = superclass and S = fine class, so the comparison is the
fixed-criterion case of the knowledge theorem: 0 <= R(D) - R(D|S) <= I(U;S),
and R(D|fine) <= R(D|coarse) <= R(D|fine) + I(U; fine | coarse).

Every rate is a bracket [dual lower bound, time-sharing upper bound] for the
declared finite problem: uniform weights on the sampled states, the listed
anchors and the listed cost matrix. It is not a bound for the population source.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp, softmax
from sklearn.linear_model import LogisticRegression
from threadpoolctl import threadpool_limits

from manifold_rd import bounds_at_distortion, conditional_ba, finite_ba, posterior_kl_bits
from semasia_cifar100 import MOBILE, REPO, RESNET, REVISION, features, metadata, shard_file

ROOT = Path(__file__).resolve().parent
RATES = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)  # Unconditional rates at which the gain is reported.


def entropy_bits(labels):
    counts = np.unique(labels, return_counts=True)[1]
    p = counts / counts.sum()
    return float(-(p * np.log2(p)).sum())


def calibrated_head(x_fit, y_fit, x_cal, y_cal, classes):
    """Multinomial head with one temperature fit on separate training images."""
    head = LogisticRegression(C=10.0, max_iter=400, tol=1e-6).fit(x_fit, y_fit)
    if not np.array_equal(head.classes_, np.arange(classes)):
        raise ValueError("the fitting subset must contain every class")

    def scores(x):
        z = head.decision_function(x)
        return np.column_stack((np.zeros(len(z)), z)) if z.ndim == 1 else z  # Binary heads return one column.

    logits = scores(x_cal)

    def loss(log_temp):
        scores = logits / np.exp(log_temp)
        return float(np.mean(logsumexp(scores, axis=1) - scores[np.arange(len(y_cal)), y_cal]))

    fitted = minimize_scalar(loss, bounds=(np.log(0.1), np.log(10)), method="bounded")
    if not fitted.success:
        raise RuntimeError("temperature fit failed")
    temp = float(np.exp(fitted.x))
    return (lambda x: softmax(scores(x) / temp, axis=1)), temp


def log_loss_bits(q, labels):
    return float(-np.log2(q[np.arange(len(labels)), labels]).mean())


def sweep(cost, contexts, betas):
    """Operating points of one finite problem; contexts=None means no shared S."""
    uniform = np.full(len(cost), 1.0 / len(cost))
    if contexts is None:
        return [finite_ba(uniform, cost, b) for b in betas]
    values, counts = np.unique(contexts, return_counts=True)
    parts = [cost[contexts == v] for v in values]
    masses = [np.full(len(c), 1.0 / len(c)) for c in parts]
    return [conditional_ba(counts / counts.sum(), masses, parts, b) for b in betas]


def gain_table(cost, knowledge, betas):
    """Brackets at the distortions where the unconditional rate equals RATES."""
    base = sweep(cost, None, betas)
    curves = {name: sweep(cost, labels, betas) for name, labels in knowledge.items()}
    rows = []
    for rate in RATES:
        if not base[0].rate < rate < base[-1].rate:
            continue
        # The unconditional achievable curve is decreasing in D: interpolate its inverse.
        d = float(np.interp(rate, [a.rate for a in base], [a.distortion for a in base]))
        row = dict(distortion=d, no_knowledge=bounds_at_distortion(base, d))
        for name, points in curves.items():
            lo, hi = bounds_at_distortion(points, d)
            row[name] = (lo, hi)
            # Outer bracket of the gain R(D) - R(D|S).
            row[f"gain_{name}"] = (row["no_knowledge"][0] - hi, row["no_knowledge"][1] - lo)
        rows.append(row)
    worst = max(a.gap for points in (base, *curves.values()) for a in points)
    return dict(rows=rows, max_lagrangian_gap_bits=worst, slopes=len(betas),
                zero_rate_distortion={"no_knowledge": float(cost.mean(0).min()), **{
                    name: float(sum(np.mean(labels == v) * cost[labels == v].mean(0).min()
                                    for v in np.unique(labels)))
                    for name, labels in knowledge.items()}})


def run(directory, download, seed, n_source, n_anchor):
    rng = np.random.default_rng(seed)
    n_fit, n_cal = 8192, 2048
    report = dict(dataset=REPO, revision=REVISION, seed=seed, fit=n_fit, calibration=n_cal,
                  source_states=n_source, anchors=n_anchor, target_unconditional_rates_bits=RATES,
                  knowledge=dict(theta="frozen backbone and calibrated fine-class head",
                                 S=["coarse_label (20)", "fine_label (100)"]), models={})
    reference, selection = {}, {}
    for model in (RESNET, MOBILE):
        arrays, labels = {}, {}
        for split, count in (("train", n_fit + n_cal), ("test", n_source)):
            shard = shard_file(directory, model, split, download)
            meta = metadata(shard, model)
            if split not in reference:
                reference[split] = meta
                selection[split] = rng.choice(len(meta), count, replace=False)
            np.testing.assert_array_equal(meta, reference[split])
            chosen = meta[selection[split]]
            arrays[split], labels[split] = features(shard, chosen[:, 0]), chosen[:, 1:]
        train, test = arrays["train"], arrays["test"]
        fine_tr, coarse_tr = labels["train"].T
        fine, coarse = labels["test"].T
        fit, cal = slice(0, n_fit), slice(n_fit, None)
        # Anchors are the first training IDs, shared by every knowledge condition.
        anchors = train[:n_anchor]

        fine_head, temp = calibrated_head(train[fit], fine_tr[fit], train[cal], fine_tr[cal], 100)
        coarse_head, _ = calibrated_head(train[fit], coarse_tr[fit], train[cal], coarse_tr[cal], 20)
        q, q_anchor = fine_head(test), fine_head(anchors)
        # H(S) - cross-entropy lower-bounds I(U;S) for the population; H(S) is the
        # value for the empirical cloud, whose distinct points determine S.
        h_coarse, h_fine = entropy_bits(coarse), entropy_bits(fine)
        info = dict(
            H_coarse_bits=h_coarse, H_fine_bits=h_fine, H_fine_given_coarse_bits=h_fine - h_coarse,
            I_U_coarse_lower_bits=h_coarse - log_loss_bits(coarse_head(test), coarse),
            I_U_fine_lower_bits=h_fine - log_loss_bits(q, fine),
            smallest_context=dict(coarse=int(np.unique(coarse, return_counts=True)[1].min()),
                                  fine=int(np.unique(fine, return_counts=True)[1].min())))

        chord = np.maximum(2.0 - 2.0 * test @ anchors.T, 0.0)  # Squared chord of unit features.
        task = posterior_kl_bits(q, q_anchor)                  # Excess fine-class log loss, bits.
        both = dict(coarse=coarse, fine=fine)
        report["models"][model] = dict(
            dimension=test.shape[1], temperature=temp, information=info,
            fine_accuracy=float(np.mean(q.argmax(1) == fine)),
            generalized_squared_chord=gain_table(chord, both, np.geomspace(1, 512, 28)),
            # The fine label is the task variable, so only the superclass is shared here.
            task_posterior_kl=gain_table(task, dict(coarse=coarse), np.geomspace(0.125, 64, 28)))
    report["id_sha256"] = {
        split: hashlib.sha256(reference[split][rows, 0].astype("<u4").tobytes()).hexdigest()
        for split, rows in selection.items()}
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="local cache of the pinned Parquet files")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--seed", type=int, default=20261116)
    parser.add_argument("--source", type=int, default=8192, help="official test images used as source states")
    parser.add_argument("--anchors", type=int, default=256)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "knowledge_cifar100.json")
    args = parser.parse_args()
    with threadpool_limits(limits=4):
        result = run(args.data, args.download, args.seed, args.source, args.anchors)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for model, entry in result["models"].items():
        print(model, json.dumps(entry["information"]))
        for name in ("generalized_squared_chord", "task_posterior_kl"):
            for row in entry[name]["rows"]:
                print(" ", name, {k: np.round(v, 4).tolist() for k, v in row.items()})
