"""Reproduce descriptor synchronization and integer precision allocation.

Run from the repository root:
  python code/synthetic/sync_allocation.py --output code/results/paper2_sync

The descriptor experiment uses actual endpoint-list and enumerative encoders,
decodes their messages, and measures bottleneck error. The allocation experiment
solves the specified source-format budget exactly; its costs are design
calculations, not measurements of a trained communication system.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
from importlib.metadata import version
import json
import math
import os
import platform
from pathlib import Path
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/semantic-manifold-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from persim import bottleneck


def ceil_log2_count(count: int) -> int:
    return (int(count) - 1).bit_length()


def list_encode(diagram: np.ndarray, precision: int, max_intervals: int,
                scale: float = 1.0) -> str:
    levels = 1 << precision
    indices = np.minimum(np.floor(diagram / scale * levels).astype(int),
                         levels - 1)
    order = np.lexsort((indices[:, 1], indices[:, 0]))
    width = ceil_log2_count(max_intervals + 1)
    message = format(len(diagram), f"0{width}b")
    for pair in indices[order]:
        message += format(int(pair[0]), f"0{precision}b")
        message += format(int(pair[1]), f"0{precision}b")
    return message


def list_decode(message: str, precision: int, max_intervals: int,
                scale: float = 1.0) -> np.ndarray:
    width = ceil_log2_count(max_intervals + 1)
    m = int(message[:width], 2)
    if len(message) != width + 2 * m * precision:
        raise ValueError("Incorrect message length")
    pairs = []
    for i in range(m):
        start = width + 2 * i * precision
        u = int(message[start:start + precision], 2)
        v = int(message[start + precision:start + 2 * precision], 2)
        if u < v:
            pairs.append(((u + .5) * scale / (1 << precision),
                          (v + .5) * scale / (1 << precision)))
    return np.asarray(pairs, dtype=float).reshape(-1, 2)


def pair_category(u: int, v: int, grid_steps: int) -> int:
    return 1 + u * grid_steps - u * (u - 1) // 2 + v - u - 1


def category_pair(category: int, grid_steps: int) -> tuple[int, int]:
    offset = category - 1
    for u in range(grid_steps):
        count = grid_steps - u
        if offset < count:
            return u, u + 1 + offset
        offset -= count
    raise ValueError("Invalid interval category")


def enumerative_encode(diagram: np.ndarray, precision: int,
                       max_intervals: int, scale: float = 1.0) -> tuple[str, int]:
    # Grid spacing is 2*epsilon, epsilon=T*2**(-s-1); include both endpoints.
    steps = 1 << precision
    indices = np.clip(np.floor(diagram / scale * steps + .5).astype(int),
                      0, steps)
    categories = [pair_category(int(u), int(v), steps)
                  for u, v in indices if u < v]
    categories += [0] * (max_intervals - len(categories))
    categories.sort()
    rank = sum(math.comb(category + j, j + 1)
               for j, category in enumerate(categories))
    n_categories = steps * (steps + 1) // 2
    possibilities = math.comb(n_categories + max_intervals, max_intervals)
    length = ceil_log2_count(possibilities)
    assert 0 <= rank < possibilities
    return format(rank, f"0{length}b"), possibilities


def enumerative_decode(message: str, precision: int, max_intervals: int,
                       scale: float = 1.0) -> np.ndarray:
    rank = int(message, 2)
    steps = 1 << precision
    n_categories = steps * (steps + 1) // 2
    upper = n_categories + max_intervals - 1
    selected = [0] * max_intervals
    for j in range(max_intervals, 0, -1):
        lo, hi = j - 1, upper
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if math.comb(mid, j) <= rank:
                lo = mid
            else:
                hi = mid - 1
        selected[j - 1] = lo
        rank -= math.comb(lo, j)
        upper = lo - 1
    assert rank == 0
    categories = [selected[j] - j for j in range(max_intervals)]
    pairs = [category_pair(c, steps) for c in categories if c > 0]
    return (np.asarray(pairs, dtype=float).reshape(-1, 2)
            * scale / steps)


def cost_and_errors(b: int, s: int, n_coordinates: int, m: int, reuse: int,
                    cx: float, ck: float) -> dict:
    header = ceil_log2_count(m + 1)
    return {"b": b, "s": s,
            "payload_error_bound": 2 * cx * 2.0 ** (-b),
            "knowledge_error_bound": ck * 2.0 ** (-s),
            "payload_bits": n_coordinates * b,
            "knowledge_bits": header + 2 * m * s,
            "amortized_bits": n_coordinates * b + (header + 2 * m * s) / reuse}


def min_precision(constant: float, budget: float) -> int:
    if budget <= 0:
        raise ValueError("Positive error budget required")
    answer = max(1, int(math.ceil(math.log2(constant / budget))))
    while constant * 2.0 ** (-answer) > budget:
        answer += 1
    return answer


def allocation(error_budget: float, n_coordinates: int, m: int, reuse: int,
               cx: float, ck: float) -> dict:
    b0 = min_precision(2 * cx, error_budget / 2)
    s0 = min_precision(ck, error_budget / 2)
    equal = cost_and_errors(b0, s0, n_coordinates, m, reuse, cx, ck)
    best = equal
    upper_b = int(math.floor(equal["amortized_bits"] / n_coordinates))
    for b in range(1, upper_b + 1):
        remainder = error_budget - 2 * cx * 2.0 ** (-b)
        if remainder <= 0:
            continue
        s = min_precision(ck, remainder)
        candidate = cost_and_errors(b, s, n_coordinates, m, reuse, cx, ck)
        if candidate["amortized_bits"] < best["amortized_bits"]:
            best = candidate
    # Independent exhaustive finite check with a provably sufficient s limit:
    # no candidate above initial objective can improve it.
    upper_s = int(math.floor(equal["amortized_bits"] * reuse / (2 * m)))
    exhaustive = float("inf")
    for b in range(1, upper_b + 1):
        for s in range(1, upper_s + 1):
            if 2 * cx * 2.0 ** (-b) + ck * 2.0 ** (-s) <= error_budget:
                exhaustive = min(exhaustive,
                                 cost_and_errors(b, s, n_coordinates, m,
                                                 reuse, cx, ck)["amortized_bits"])
    assert math.isclose(exhaustive, best["amortized_bits"])
    ax, ak = n_coordinates, 2 * m / reuse
    ex, ek = ax * error_budget / (ax + ak), ak * error_budget / (ax + ak)
    bstar, sstar = math.log2(2 * cx / ex), math.log2(ck / ek)
    # All supplied configurations lie in the interior regime.
    assert bstar >= 1 and sstar >= 1
    relaxed = ax * bstar + ak * sstar + ceil_log2_count(m + 1) / reuse
    assert best["amortized_bits"] + 1e-9 >= relaxed
    return {"error_budget": error_budget, "n_coordinates": n_coordinates,
            "knowledge_intervals": m, "reuse": reuse,
            "equal_error_b": equal["b"], "equal_error_s": equal["s"],
            "equal_error_bits": equal["amortized_bits"],
            "optimal_b": best["b"], "optimal_s": best["s"],
            "optimal_bits": best["amortized_bits"],
            "continuous_b": bstar, "continuous_s": sstar,
            "continuous_bits": relaxed,
            "saving_fraction": 1 - best["amortized_bits"] / equal["amortized_bits"],
            "achieved_budget": best["payload_error_bound"] + best["knowledge_error_bound"]}


def write_csv(path: Path, records: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("code/results/paper2_sync"))
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--trials", type=int, default=20)
    args = parser.parse_args()
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    records, summaries = [], []
    for m in (4, 8, 16):
        # A separated long-lifetime family plus variable cardinality tests
        # exercise both diagram geometry and count coding.
        diagrams = []
        for trial in range(args.trials):
            cardinality = m if trial < args.trials - 2 else (0 if trial == args.trials - 1 else 1)
            diagrams.append(np.column_stack((rng.uniform(0, .4, cardinality),
                                            rng.uniform(.6, 1, cardinality))))
        for precision in range(2, 11):
            epsilon = 2.0 ** (-precision - 1)
            for trial, diagram in enumerate(diagrams):
                list_message = list_encode(diagram, precision, m)
                list_rec = list_decode(list_message, precision, m)
                enum_message, possibilities = enumerative_encode(diagram, precision, m)
                enum_rec = enumerative_decode(enum_message, precision, m)
                list_error = float(bottleneck(diagram, list_rec)) if len(diagram) else 0.0
                enum_error = float(bottleneck(diagram, enum_rec)) if len(diagram) else 0.0
                assert list_error <= epsilon + 1e-12
                assert enum_error <= epsilon + 1e-12
                expected_list_bits = ceil_log2_count(m + 1) + 2 * len(diagram) * precision
                assert len(list_message) == expected_list_bits
                assert len(enum_message) == ceil_log2_count(possibilities)
                records.append({"max_intervals": m, "cardinality": len(diagram),
                                "precision": precision, "trial": trial,
                                "epsilon": epsilon,
                                "list_bits": len(list_message),
                                "enumerative_bits": len(enum_message),
                                "list_bottleneck": list_error,
                                "enumerative_bottleneck": enum_error})
            full = [row for row in records if row["max_intervals"] == m
                    and row["precision"] == precision and row["cardinality"] == m]
            summaries.append({"max_intervals": m, "precision": precision,
                              "epsilon": epsilon,
                              "list_bits": full[0]["list_bits"],
                              "enumerative_bits": full[0]["enumerative_bits"],
                              "mean_list_bottleneck": float(np.mean([r["list_bottleneck"] for r in full])),
                              "max_list_bottleneck": max(r["list_bottleneck"] for r in full),
                              "max_enumerative_bottleneck": max(r["enumerative_bottleneck"] for r in full)})
    cx, ck, ncoord, m = math.sqrt(2) * 2.5, .5, 256, 32
    budgets = sorted(set(float(x) for x in np.geomspace(.02, .25, 25)) | {.111})
    allocations = [allocation(budget, ncoord, m, reuse, cx, ck)
                   for reuse in (1, 4, 16, 64) for budget in budgets]
    examples = [r for r in allocations if r["error_budget"] == .111]
    write_csv(args.output / "coding_trials.csv", records)
    write_csv(args.output / "synchronization.csv", summaries)
    write_csv(args.output / "allocations.csv", allocations)
    write_csv(args.output / "allocation_examples.csv", examples)
    metadata = {"seed": args.seed, "trials_per_cardinality": args.trials,
                "scale_range": 1.0, "n_coordinates": ncoord,
                "coordinate_quantizer_half_range": 2.5,
                "allocation_knowledge_intervals": m,
                "allocation_residual_error_budget_example": .111,
                "allocation_link_model": "digital delivery conditioned on successful source-message decoding",
                "python": platform.python_version(), "numpy": np.__version__,
                "persim": version("persim"), "platform": platform.platform(),
                "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "command": "python " + " ".join(sys.argv),
                "elapsed_seconds": time.monotonic() - started}
    report = {"metadata": metadata,
              "checks": {"coded_reconstructions": 2 * len(records),
                         "endpoint_bound_violations": 0,
                         "encoding_length_violations": 0,
                         "allocation_configurations": len(allocations),
                         "integer_optimum_violations": 0},
              "coding_summaries": summaries, "allocation_examples": examples,
              "median_saving_fraction_all_budgets": {
                  str(reuse): float(np.median([r["saving_fraction"] for r in allocations
                                               if r["reuse"] == reuse]))
                  for reuse in (1, 4, 16, 64)}}
    (args.output / "results.json").write_text(json.dumps(report, indent=2) + "\n",
                                              encoding="utf-8")
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.0), layout="constrained")
    subset = [r for r in summaries if r["max_intervals"] == 16]
    ss = [r["precision"] for r in subset]
    axes[0].semilogy(ss, [r["epsilon"] for r in subset], "k--", label="Worst-case bound")
    axes[0].semilogy(ss, [r["max_list_bottleneck"] for r in subset], "o-", color="#126682",
                     label="List: maximum measured")
    axes[0].semilogy(ss, [r["max_enumerative_bottleneck"] for r in subset], "s:", color="#b86b26",
                     label="Enumerative: maximum measured")
    axes[0].set(xlabel="Endpoint precision s (bits)", ylabel="Bottleneck error",
                title="(a) Descriptor reconstruction")
    axes[0].legend(fontsize=7)
    axes[1].plot(ss, [r["list_bits"] for r in subset], "o-", color="#126682", label="Endpoint list")
    axes[1].plot(ss, [r["enumerative_bits"] for r in subset], "s--", color="#b86b26", label="Enumerative code")
    axes[1].set(xlabel="Equivalent endpoint precision (bits)", ylabel="Source length (bits)",
                title="(b) Unordered descriptor coding")
    axes[1].legend(fontsize=7)
    hh = [r["reuse"] for r in examples]
    axes[2].plot(hh, [r["equal_error_bits"] for r in examples], "o-", color="#126682",
                 label="Equal error split")
    axes[2].plot(hh, [r["optimal_bits"] for r in examples], "s--", color="#b86b26",
                 label="Optimal integer allocation")
    axes[2].plot(hh, [r["continuous_bits"] for r in examples], "^:", color="#555555",
                 label="Continuous relaxation")
    axes[2].set_xscale("log", base=4)
    axes[2].set_xticks(hh, [str(h) for h in hh])
    axes[2].set(xlabel="Knowledge reuse H (blocks)", ylabel="Amortized source bits / block",
                title="(c) Joint precision allocation")
    axes[2].legend(fontsize=7)
    for axis in axes:
        axis.grid(alpha=.18)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(args.output / f"synchronization_allocation.{extension}", dpi=220)
    plt.close(fig)
    print(json.dumps({"checks": report["checks"], "examples": examples,
                      "median_savings": report["median_saving_fraction_all_budgets"]}, indent=2))


if __name__ == "__main__":
    main()

