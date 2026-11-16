import math, sys, unittest
from pathlib import Path

d = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(d), str(d / "experiments")]
import numpy as np
import torch
from manifold_rd import finite_ba, conditional_ba
from semasia_rd import solve
from reference_graph import fit_graph, attached
from link import transmit


class Certificates(unittest.TestCase):
    def test_binary_against_closed_form(self):
        c = torch.tensor([[0.0, 1.0], [1.0, 0.0]], dtype=torch.float64)
        beta = 2.5
        _, v = solve(c, beta, 80, 1e-10)
        d = 1 / (1 + 2**beta)
        opt = 1 + d * math.log2(d) + (1 - d) * math.log2(1 - d) + beta * d
        self.assertLessEqual(v["lower"], opt + 1e-10)
        self.assertGreaterEqual(v["upper"], opt - 1e-10)
        self.assertAlmostEqual(v["distortion"], d, places=9)
        ck = conditional_ba(
            [0.5, 0.5], [[0.5, 0.5], [1.0, 0.0]], [c.numpy(), c.numpy()], beta
        )
        self.assertAlmostEqual(
            ck.rate, 0.5 * finite_ba([0.5, 0.5], c.numpy(), beta).rate, places=7
        )

    def test_graph_arc_and_heldout_attachment(self):
        th = torch.arange(12, dtype=torch.float64) * math.pi / 6
        z = torch.stack((th.cos(), th.sin()), 1)
        graph, c = fit_graph(z, 1, 2)
        self.assertAlmostEqual(
            float(c[6, 0].sqrt()), 12 * math.sin(math.pi / 12), places=8
        )
        h = torch.tensor(
            [[math.cos(math.pi / 12), math.sin(math.pi / 12)]], dtype=z.dtype
        )
        self.assertAlmostEqual(
            float(attached(h, graph)[0, 0].sqrt()), 2 * math.sin(math.pi / 24), places=8
        )
        with self.assertRaisesRegex(ValueError, "components"):
            fit_graph(torch.tensor([[0.0], [1.0], [10.0], [11.0]]), 1, 1)

    def test_operational_bpsk(self):
        gen = torch.Generator().manual_seed(41)
        idx = torch.arange(60000) % 16
        decoded, wrong = transmit(idx, 4, "awgn", 0, gen)
        self.assertAlmostEqual(wrong / (4 * len(idx)), math.erfc(1) / 2, delta=0.003)
        np.testing.assert_array_equal(transmit(idx, 4, "clean", None, gen)[0], idx)


if __name__ == "__main__":
    unittest.main()
