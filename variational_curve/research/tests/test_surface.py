import sys, unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import torch
from semasia_curve import Surface
from link import information, LN2


class Constraints(unittest.TestCase):
    def test_cannot_collapse_or_expand_secants(self):
        torch.manual_seed(2)
        net = Surface(11, 5, 4, 0.6, "joint")
        with torch.no_grad():
            net.a.mul_(12)
            net.b.normal_(0, 13)
        x, y = torch.randn(50, 11), torch.randn(50, 11)
        ratio = (net(x) - net(y)).norm(dim=1) / (x - y).norm(dim=1)
        self.assertTrue(torch.all((ratio >= 0.4) & (ratio <= 1.6)))

    def test_covariance_gradient_and_free_marginal(self):
        torch.manual_seed(3)
        logits = torch.randn(7, 4, dtype=torch.float64, requires_grad=True)
        W = logits.softmax(1)
        c = torch.rand_like(W)
        L = information(W) + (W * c).sum(1).mean()
        (g,) = torch.autograd.grad(L, logits)
        u = (W.log() - W.mean(0).log()) / LN2 + c
        torch.testing.assert_close(g, W * (u - (W * u).sum(1, keepdim=True)) / len(W))
        obj = (W * (W.log() - torch.zeros(4, dtype=W.dtype).log_softmax(0))).sum(
            1
        ).mean() / LN2
        self.assertGreaterEqual(float(obj.detach()), float(information(W).detach()))


if __name__ == "__main__":
    unittest.main()
