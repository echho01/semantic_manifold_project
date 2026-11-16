"""Regression checks for mathematical and experiment-pipeline invariants."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from threadpoolctl import threadpool_limits
CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE / "embeddings"))
sys.path.insert(0, str(CODE / "synthetic"))
import geometry
import manifold_gain
import rd
import task_mismatch
from run_paper1 import sphere_rd, projection_sphere, surface_geometry, surface_ba


class GeometryTests(unittest.TestCase):
    def test_high_ambient_curvature_and_translation(self):
        rng = np.random.default_rng(4)
        x = rng.normal(size=(4000, 3))
        x = x / np.linalg.norm(x, axis=1, keepdims=True) * 3
        q, _ = np.linalg.qr(rng.normal(size=(80, 3)))
        x = x @ q.T
        original = geometry.local_curvature(x, 2, k_nn=60, n_anchor=40)
        translated = geometry.local_curvature(x + 3.2, 2, k_nn=60, n_anchor=40)
        self.assertLess(abs(original["scal_mean"] / (2 / 9) - 1), .06)
        self.assertAlmostEqual(original["scal_mean"], translated["scal_mean"], places=10)
        self.assertAlmostEqual(original["h2_mean"] - original["ii2_mean"], original["scal_mean"], places=12)

    def test_affine_decoder_centres_and_tangent_dimension(self):
        rng = np.random.default_rng(9)
        q, _ = np.linalg.qr(rng.normal(size=(32, 2)))
        reference = rng.normal(size=(2000, 2)) @ q.T + 5
        source = rng.normal(size=(400, 2)) @ q.T + 5
        noise = rng.normal(scale=.01, size=source.shape)
        decoded = manifold_gain.tangent_project(source + noise, reference, k=24, pca_dim=2)
        truth = source + (noise @ q) @ q.T
        np.testing.assert_allclose(decoded, truth, atol=1e-11)

    def test_twonn_rejects_discrete_prompts(self):
        x = np.repeat(np.eye(10), 50, axis=0)
        with self.assertRaisesRegex(ValueError, "distinct"):
            geometry.twonn_dimension(x)


class RateTests(unittest.TestCase):
    def test_continuous_sphere_coefficients(self):
        for d, criterion, prediction in ((2, "chordal", 0.), (3, "geodesic", 1/3), (1, "chordal", -.125)):
            rows = [sphere_rd(d, mu, criterion) for mu in (1000, 2000)]
            fitted = np.polyfit([r["D"] for r in rows], [r["residual_over_D"] for r in rows], 1)[1]
            self.assertAlmostEqual(fitted, prediction, places=6)

    def test_projection_negative_second_order(self):
        variance = .001
        coefficient = (projection_sphere(3, variance) - 3*variance)/variance**2
        self.assertLess(abs(coefficient+1), .01)

    def test_positive_ellipsoid_area_and_ba_gap(self):
        geom = surface_geometry("ellipsoid", n=160)
        self.assertTrue(np.all(geom["mass"] > 0))
        self.assertAlmostEqual(geom["mass"].sum(), 1., places=12)
        torus = surface_geometry("torus", n=256)
        result = surface_ba(torus, mu=100, tolerance=1e-10)
        self.assertLessEqual(result["lagrangian_gap_nats"], 1e-10)
        self.assertLess(abs(result["residual_over_D"] - torus["coefficient"]), .001)

    def test_heldout_quantizer_and_unsupported_arc(self):
        x = np.random.default_rng(10).normal(size=(400, 4))
        curve = rd.quantize_curve(x, [2, 4, 8, 1000], seed=1)
        self.assertEqual(curve["train_count"], 320)
        self.assertEqual(curve["test_count"], 80)
        self.assertEqual(curve["sizes"], [2, 4, 8])
        self.assertGreater(curve["D"][0], curve["D"][-1])
        json.dumps(curve, allow_nan=False)
        with self.assertRaisesRegex(ValueError, "geodesic"):
            rd.quantize_curve(x, [2], geodesic=True, tau=1.)


class TaskTests(unittest.TestCase):
    def test_probe_dimension_and_same_edge_metric(self):
        rng = np.random.default_rng(12)
        x = rng.normal(size=(500, 8))
        y = np.arange(500) % 4
        probe = task_mismatch.linear_probe(x, y, 2)
        self.assertEqual(probe.shape, (2, 8))
        self.assertEqual(np.linalg.matrix_rank(probe), 2)
        eps = task_mismatch.estimate_epsilon(x, 2*x, k=8)
        self.assertAlmostEqual(eps["median"], 3., places=12)
        self.assertAlmostEqual(eps["max_observed"], 3., places=12)

    def test_transfer_preserves_the_source_index(self):
        uv = np.random.default_rng(14).normal(size=(2000, 2))
        a = uv * [100, 1]
        b = uv * [1, 100]
        result = task_mismatch.mismatch_rd(a, b, [8], seed=4)
        self.assertGreater(result["transfer_over_matched_B"][0], 5.)


class DriverTests(unittest.TestCase):
    def test_cpu_stages_resume_and_outage_only(self):
        rng = np.random.default_rng(31)
        z = rng.normal(size=(600, 3))
        z /= np.linalg.norm(z, axis=1, keepdims=True)
        q, _ = np.linalg.qr(rng.normal(size=(12, 3)))
        x = z @ q.T
        labels = (z[:,0] > 0).astype(int)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for model, embeddings in (("a", x), ("b", 1.2*x)):
                run = root / "runs" / f"{model}__toy"
                run.mkdir(parents=True)
                np.savez(run / "arrays.npz", image=embeddings, labels=labels, item_ids=np.arange(len(x)))
            cfg = {"seed": 0, "paths": {"root": "runs", "hf_cache": "hf", "torch_cache": "torch"},
                   "models": {"a": {}, "b": {}}, "datasets": {"toy": {}},
                   "geometry": {"twonn": {"max_points": 600, "frac": .9},
                                "curvature": {"pca_dim": 2, "k_nn": 32, "n_anchor": 10}},
                   "rd": {"codebook_sizes": [2,4,8], "criteria": ["chordal"], "kmeans_iters": 5},
                   "gain_outage": {"noise_levels": [.01,.1], "n_points": 30, "n_neighbors": 16},
                   "task": {"probe_dims": [1,2], "codebook_sizes": [2,4]},
                   "mismatch": {"model_a": "a", "model_b": "b", "codebook_sizes": [2,4]}}
            config = root / "config.json"
            config.write_text(json.dumps(cfg))
            command = [sys.executable, str(CODE / "embeddings" / "run_all.py"), "--config", str(config), "--model", "a", "--dataset", "toy"]
            subprocess.run(command + ["--stages", "geometry,rd,gain,outage,task,mismatch"], check=True, capture_output=True, text=True)
            subprocess.run(command + ["--stages", "outage"], check=True, capture_output=True, text=True)
            report = json.loads((root / "runs/a__toy/report.json").read_text())
            self.assertIn("rd", report)
            self.assertEqual(len(report["runs"]), 2)
            self.assertEqual(report["gain"]["pca_dim"], 2)
            self.assertIn("curve", report["mismatch"])


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        unittest.main()
