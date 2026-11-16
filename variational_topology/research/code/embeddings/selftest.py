"""Self-test of the geometry estimators on manifolds with known invariants.

Run:  python selftest.py
"""
import numpy as np
import geometry


def sphere(n, d, D, rho, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, d + 1))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    X = np.zeros((n, D), np.float32)
    X[:, :d + 1] = rho * v
    return X


def main():
    print("== TwoNN on flat clouds (expected: the ambient dimension) ==")
    rng = np.random.default_rng(1)
    for D in [2, 5, 20]:
        X = rng.normal(size=(20000, D)).astype(np.float32)
        r = geometry.twonn_dimension(X, max_points=20000, frac=0.9)
        print(f"  R^{D}: d_hat={r['dim']:.3f} residual={r['residual']:.3f}")

    print("== TwoNN on a 2-plane in R^50 (expected 2) ==")
    n = 20000
    Y = rng.normal(size=(n, 2)) @ rng.normal(size=(2, 50))
    print("  ", geometry.twonn_dimension(Y.astype(np.float32), max_points=n))

    print("== TwoNN and curvature on S^2(rho=3) in R^10 ==")
    X = sphere(8000, 2, 10, 3.0)
    print("  TwoNN:", geometry.twonn_dimension(X, max_points=8000))
    c = geometry.local_curvature(X, pca_dim=2, k_nn=80, n_anchor=800, seed=0)
    print(f"  Scal mean = {c['scal_mean']:.4f}   expected 2/rho^2 = {2/9:.4f}")
    print(f"  reach med = {c['reach_median']:.3f}  expected rho = 3")
    print(f"  |II|^2 mean = {c['ii2_mean']:.4f}  expected d/rho^2 = {2/9:.4f}")


if __name__ == "__main__":
    main()
