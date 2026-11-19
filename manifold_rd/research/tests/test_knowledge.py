"""Fixed-criterion knowledge gain on a synthetic clustered cloud (NumPy/SciPy/scikit-learn only)."""
import sys, unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from knowledge_cifar100 import entropy_bits, gain_table


class KnowledgeGain(unittest.TestCase):
    def test_side_information_bounds(self):
        rng = np.random.default_rng(7)
        n, anchors, d = 600, 32, 8
        coarse = rng.integers(0, 4, n)
        fine = coarse * 3 + rng.integers(0, 3, n)
        centers = rng.normal(size=(12, d))
        unit = lambda v: v / np.linalg.norm(v, axis=1, keepdims=True)
        x = unit(centers[fine] + 0.4 * rng.normal(size=(n, d)))
        a = unit(centers[rng.integers(0, 12, anchors)] + 0.4 * rng.normal(size=(anchors, d)))
        cost = np.maximum(2 - 2 * x @ a.T, 0)
        table = gain_table(cost, dict(coarse=coarse, fine=fine), np.geomspace(1, 256, 20))
        self.assertGreaterEqual(len(table["rows"]), 3)
        slack = 1e-6
        for row in table["rows"]:
            for name, labels in (("coarse", coarse), ("fine", fine)):
                low, high = row[f"gain_{name}"]
                # 0 <= R(D) - R(D|S) <= I(U;S) = H(S) for a cloud of distinct points.
                self.assertGreaterEqual(high, -slack)
                self.assertLessEqual(low, entropy_bits(labels) + slack)
            # Nested knowledge: R(D|fine) <= R(D|coarse) <= R(D|fine) + H(fine|coarse).
            self.assertLessEqual(row["fine"][0], row["coarse"][1] + slack)
            self.assertLessEqual(row["coarse"][0],
                                 row["fine"][1] + entropy_bits(fine) - entropy_bits(coarse) + slack)
        zero = table["zero_rate_distortion"]
        self.assertLessEqual(zero["fine"], zero["coarse"] + slack)
        self.assertLessEqual(zero["coarse"], zero["no_knowledge"] + slack)


if __name__ == "__main__":
    unittest.main()
