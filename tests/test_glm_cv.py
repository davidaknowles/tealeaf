"""Checks for scalable GLM count-fold cross-validation."""

import unittest

import numpy as np
import scipy.sparse as sp

try:
    from tealeaf.sc import glm_cv, glm_solvers
    glm_solvers._torch()
except ImportError:  # pragma: no cover
    TORCH_AVAILABLE = False
else:
    TORCH_AVAILABLE = True


@unittest.skipUnless(TORCH_AVAILABLE, "Torch optional dependency is unavailable")
class GLMCVTest(unittest.TestCase):
    def setUp(self):
        self.counts = sp.csr_matrix(
            np.array(
                [[20, 3, 2], [2, 18, 4], [8, 9, 3], [5, 4, 11]],
                dtype=np.int64,
            )
        )
        self.compatibility = sp.csr_matrix(
            np.array(
                [[1.0, 0.1], [0.1, 1.0], [0.7, 0.7]], dtype=np.float32
            )
        )

    def test_count_folds_partition_input(self):
        folds = glm_cv.split_count_folds(self.counts, n_folds=3, seed=4)
        self.assertEqual(len(folds), 3)
        reconstructed = sum(folds, start=sp.csr_matrix(self.counts.shape))
        np.testing.assert_array_equal(reconstructed.toarray(), self.counts.toarray())
        repeated = glm_cv.split_count_folds(self.counts, n_folds=3, seed=4)
        for observed, expected in zip(folds, repeated):
            np.testing.assert_array_equal(observed.toarray(), expected.toarray())

    def test_hyperparameter_scales_are_positive(self):
        for method in ("admm_factorized", "frank_wolfe_penalized"):
            scale = glm_cv.hyperparameter_scale(
                self.counts,
                self.compatibility,
                method,
                device="cpu",
                batch_cells=2,
                power_iter=5,
            )
            self.assertTrue(np.isfinite(scale))
            self.assertGreater(scale, 0)

    def test_only_open_grid_boundaries_require_expansion(self):
        self.assertFalse(
            glm_cv._best_on_open_boundary("admm_factorized", [0, 0.1, 1], 1)
        )
        self.assertFalse(
            glm_cv._best_on_open_boundary("frank_wolfe_penalized", [0, 1, 4], 0)
        )
        self.assertTrue(
            glm_cv._best_on_open_boundary("frank_wolfe_penalized", [0, 1, 4], 4)
        )

    def test_cross_validation_reports_best_multiplier(self):
        report = glm_cv.cross_validate_glm(
            self.counts,
            self.compatibility,
            "frank_wolfe_penalized",
            [0.5, 2.0],
            n_folds=2,
            device="cpu",
            batch_cells=2,
            power_iter=3,
            fit_kwargs={
                "rank": 2,
                "max_atoms": 2,
                "max_iter": 2,
                "min_iter": 2,
                "power_iter": 3,
            },
        )
        self.assertIn(report["best_multiplier"], [0.5, 2.0])
        self.assertIsInstance(report["best_on_boundary"], bool)
        self.assertEqual(len(report["fold_results"]), 4)


if __name__ == "__main__":
    unittest.main()
