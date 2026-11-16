"""CPU self-test of the Paper 2 pipeline on synthetic 'models' with a known
answer.  No torch needed.

Two fake encoders of the same 3000 items (5 classes x angle):
  model_a  each class on its own circle (radius 0.6) around a class center
           in R^64: 5 concepts, 5 loops
  model_b  the same items, but classes 3 and 4 share one circle (a receiver
           KB that merges two concepts): a 'language' mismatch
Expected: floor(a,b) > 0 and d_topo(K) -> floor within eps_R(K); the
reparametrization bounds hold; the raw receiver's spurious onset sits at
sqrt2*theta; the debiased count stays at 5 well beyond it; the KB-aware
receiver is the most robust.

Run:  python selftest_topology.py
"""
import json
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

# macOS Accelerate emits spurious matmul warnings; ripser warns on wide inputs
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

HERE = Path(__file__).parent


def fake_models(root: Path, n: int = 3000, seed: int = 0):
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 5, n)
    th = rng.uniform(0, 2 * np.pi, n)
    C = rng.normal(size=(5, 64))
    C *= 4.0 / np.linalg.norm(C, axis=1, keepdims=True)
    F = np.linalg.qr(rng.normal(size=(64, 2)))[0].T            # one shared circle plane
    A = C[y] + 0.6 * np.c_[np.cos(th), np.sin(th)] @ F
    yb = np.where(y == 4, 3, y)                                # classes 3,4 merged in b
    thb = np.where(y >= 3, th / 2 + np.pi * (y == 4), th)      # and share one circle
    B = C[yb] + 0.6 * np.c_[np.cos(thb), np.sin(thb)] @ F
    for name, X in (("model_a", A), ("model_b", B)):
        d = root / f"{name}__synth"
        d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / "arrays.npz", image=X.astype(np.float32), labels=y, item_ids=np.arange(n))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "runs"
        fake_models(root)
        cfg = {"seed": 0, "paths": {"root": str(root)},
               "topology": {"n_items": 240, "maxdim": 1,
                            "dictionary_sizes": [32, 64, 128, 240],
                            "reparam_points": 120, "noise_levels": [0.25, 1.0],
                            "cliff_noise_levels": [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2],
                            "cliff_per_class": 60, "reeb_bins": 20}}
        cpath = Path(tmp) / "cfg.json"
        cpath.write_text(json.dumps(cfg))
        subprocess.run([sys.executable, str(HERE / "run_topology.py"), "--config", str(cpath),
                        "--dataset", "synth", "--models", "model_a,model_b"],
                       check=True, stdout=subprocess.DEVNULL)
        rep = json.loads((root / "topology__synth" / "report.json").read_text())

    ok = True
    fl = rep["floor"]["floor"]["model_a"]["model_b"]
    print(f"floor d_B(a, b) = {fl:.4f}  (> 0: language mismatch)")
    ok &= fl > 0.05
    dc = rep["dictionary"]["model_a->model_b"]
    print(f"dictionary a->b: d_topo = {np.round(dc['d_topo'], 4).tolist()}  floor {dc['floor']:.4f}")
    print(f"                 eps_R  = {np.round(dc['eps_R'], 4).tolist()}  "
          f"within bound: {all(dc['within_bound'])}")
    ok &= all(dc["within_bound"])
    rp = rep["reparam"]["model_a"]
    print(f"reparam bounds hold on model_a: {rp['all_within']}")
    ok &= rp["all_within"]
    cl = rep["cliff"]["model_a"]
    print(f"cliff model_a: theta={cl['theta']:.3f}, clean count {cl['clean_count']}, "
          f"pred spurious sigma {cl['pred_spurious_sigma']:.3f}")
    for r in cl["rows"]:
        print(f"   sigma={r['sigma']:<5} raw={r['raw_count']:<4} deb={r['debiased_count']:<3} "
              f"KB-aware acc={r['kb_aware_acc']:.3f}")
    first_spur = next((r["sigma"] for r in cl["rows"] if r["raw_count"] > cl["clean_count"]), None)
    ok &= cl["clean_count"] == 5 and first_spur is not None
    rb = rep["reeb"]["model_a"]
    print(f"reeb model_a: H(E|f) = {rb['H_E_given_f_bits']:.3f} bits "
          f"(<= log2 max comps {rb['log2_max_components']:.3f}); mapper {rb['mapper']}")
    ok &= rb["H_E_given_f_bits"] <= rb["log2_max_components"] + 1e-9
    print("ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
