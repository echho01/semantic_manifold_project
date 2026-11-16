import sys, unittest, itertools
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import numpy as np
import torch
from scipy.ndimage import label
from structure import merge_values, cubical, w2
from topology_jscc import Codec


class Filtration(unittest.TestCase):
    def test_threshold_oracle_and_risk(self):
        rng = np.random.default_rng(312)
        f = rng.uniform(size=(4, 5))
        g = rng.uniform(size=(4, 5))
        pairs = np.array(list(itertools.product(range(20), repeat=2)))
        ts = np.unique(np.r_[0, 1, f.ravel(), g.ravel()])
        mu = merge_values(torch.tensor(f), pairs).numpy()
        nu = merge_values(torch.tensor(g), pairs).numpy()
        integral = np.zeros_like(mu)
        for lo, hi in zip(ts[:-1], ts[1:]):
            t = (lo + hi) / 2
            active = []
            for a in (f, g):
                cc = label(a <= t)[0].ravel()
                u, v = cc[pairs[:, 0]], cc[pairs[:, 1]]
                active.append((u > 0) & (v > 0) & (u == v))
            np.testing.assert_array_equal(active[0], mu <= t)
            integral += (hi - lo) * (active[0] ^ active[1])
        np.testing.assert_allclose(integral, np.abs(mu - nu), atol=1e-12)

    def test_hole_and_diagonal_gradient(self):
        v = torch.zeros((5, 5), dtype=torch.float64)
        v[2, 2] = 1
        v.requires_grad_()
        _, h1 = cubical(v)
        self.assertEqual(len(h1), 1)
        d = w2(h1, v.new_empty((0, 2)))
        self.assertAlmostEqual(float(d.detach()), 0.5)
        d.backward()
        self.assertAlmostEqual(float(v.grad[2, 2]), 0.5)

    def test_energy_for_learned_unequal_groups(self):
        net = Codec(16, 0.2, 4, 3, "learned")
        net.log_power.data.copy_(torch.tensor([-2.0, 0.0, 1.0]))
        _, z, _ = net(torch.rand(3, 1, 16, 16), "rayleigh", 0)
        e = z.reshape(3, -1, 2).square().sum(2).mean(1)
        torch.testing.assert_close(e, torch.ones_like(e))


if __name__ == "__main__":
    unittest.main()
