import unittest

import numpy as np
import pandas as pd

from extra_scripts.select_path_smoothing_from_pvalues import (
    beta_uniform_mixture,
    bh_count,
)


class TestPathSmoothingSelection(unittest.TestCase):
    def test_beta_uniform_mixture_favors_small_p_values(self):
        rng = np.random.default_rng(3)
        null = pd.DataFrame({
            "gene_id": np.arange(1000),
            "p_value": rng.uniform(size=1000),
        })
        signal = null.copy()
        signal.loc[:199, "p_value"] = rng.beta(0.2, 1.0, size=200)
        self.assertGreater(
            beta_uniform_mixture(signal)["mean_log_evidence"],
            beta_uniform_mixture(null)["mean_log_evidence"],
        )

    def test_gene_equal_weighting_ignores_test_multiplicity(self):
        base = pd.DataFrame({
            "gene_id": ["a", "b"],
            "p_value": [0.01, 0.8],
        })
        repeated = pd.concat([base.iloc[[0]]] * 20 + [base.iloc[[1]]])
        self.assertAlmostEqual(
            beta_uniform_mixture(base)["mean_log_evidence"],
            beta_uniform_mixture(repeated)["mean_log_evidence"],
            places=8,
        )

    def test_bh_count(self):
        self.assertEqual(bh_count([0.001, 0.01, 0.8]), 2)
        self.assertEqual(bh_count([0.1, 0.2]), 0)


if __name__ == "__main__":
    unittest.main()
