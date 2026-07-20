"""Checks for scalable GLM count-fold cross-validation."""

import unittest
from unittest import mock

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
        self.assertEqual(
            glm_cv._expanded_candidate(
                "frank_wolfe_penalized", [0, 1, 4], "upper", 4
            ),
            16,
        )
        self.assertEqual(
            glm_cv._expanded_candidate(
                "admm_factorized", [0, 0.1], "upper", 100
            ),
            1,
        )

    def test_adaptive_grid_only_fits_new_boundary_candidate(self):
        initial = {
            "method": "frank_wolfe_penalized",
            "n_folds": 2,
            "multipliers": [1.0, 2.0],
            "fold_scales": [1.0, 1.0],
            "mean_validation_loss": {1.0: 0.2, 2.0: 0.1},
            "validation_standard_error": {1.0: 0.01, 2.0: 0.01},
            "candidate_converged": {1.0: True, 2.0: True},
            "candidate_nondegenerate": {1.0: True, 2.0: True},
            "mean_profile_relative_variance": {1.0: 0.1, 2.0: 0.2},
            "minimum_profile_active_fraction": {1.0: 1.0, 2.0: 1.0},
            "best_multiplier": 2.0,
            "best_on_boundary": True,
            "fold_results": [{"multiplier": 1.0}, {"multiplier": 2.0}],
        }
        extension = {
            "method": "frank_wolfe_penalized",
            "n_folds": 2,
            "multipliers": [8.0],
            "fold_scales": [1.0, 1.0],
            "mean_validation_loss": {8.0: 0.3},
            "validation_standard_error": {8.0: 0.01},
            "candidate_converged": {8.0: True},
            "candidate_nondegenerate": {8.0: True},
            "mean_profile_relative_variance": {8.0: 0.3},
            "minimum_profile_active_fraction": {8.0: 1.0},
            "best_multiplier": 8.0,
            "best_on_boundary": False,
            "fold_results": [{"multiplier": 8.0}],
        }
        with mock.patch.object(
            glm_cv, "cross_validate_glm", side_effect=[initial, extension]
        ) as mocked:
            report = glm_cv.cross_validate_glm_adaptive_grid(
                self.counts,
                self.compatibility,
                "frank_wolfe_penalized",
                [1, 2],
                max_grid_expansions=3,
                grid_expansion_factor=4,
            )
        self.assertEqual(mocked.call_count, 2)
        self.assertEqual(mocked.call_args_list[1].args[3], [8.0])
        self.assertEqual(report["best_multiplier"], 2.0)
        self.assertFalse(report["best_on_boundary"])
        self.assertFalse(report["grid_exhausted"])
        self.assertEqual(report["grid_expansions"], 1)

    def test_one_se_selects_most_regularized_converged_candidate(self):
        report = {
            "multipliers": [1.0, 4.0, 16.0],
            "mean_validation_loss": {1.0: 0.12, 4.0: 0.101, 16.0: 0.1},
            "validation_standard_error": {1.0: 0.01, 4.0: 0.005, 16.0: 0.01},
            "candidate_converged": {1.0: True, 4.0: True, 16.0: True},
            "candidate_nondegenerate": {1.0: True, 4.0: True, 16.0: True},
        }
        glm_cv._apply_selection_rule(
            report, "frank_wolfe_penalized", "one_standard_error", True
        )
        self.assertEqual(report["minimum_loss_multiplier"], 16.0)
        self.assertEqual(report["best_multiplier"], 4.0)
        report["candidate_converged"][4.0] = False
        glm_cv._apply_selection_rule(
            report, "frank_wolfe_penalized", "one_standard_error", True
        )
        self.assertEqual(report["best_multiplier"], 16.0)
        report["candidate_converged"][4.0] = True
        report["candidate_nondegenerate"][4.0] = False
        glm_cv._apply_selection_rule(
            report,
            "frank_wolfe_penalized",
            "one_standard_error",
            True,
            True,
        )
        self.assertEqual(report["best_multiplier"], 16.0)

    def test_variance_retention_avoids_over_regularized_profile(self):
        report = {
            "multipliers": [1.0, 4.0, 16.0, 64.0],
            "mean_validation_loss": {
                1.0: 0.11, 4.0: 0.105, 16.0: 0.101, 64.0: 0.1,
            },
            "validation_standard_error": {
                1.0: 0.02, 4.0: 0.02, 16.0: 0.02, 64.0: 0.02,
            },
            "candidate_converged": {
                1.0: True, 4.0: True, 16.0: True, 64.0: True,
            },
            "candidate_nondegenerate": {
                1.0: True, 4.0: True, 16.0: True, 64.0: True,
            },
            "mean_profile_relative_variance": {
                1.0: 0.2, 4.0: 0.91, 16.0: 1.0, 64.0: 0.98,
            },
        }
        glm_cv._apply_selection_rule(
            report,
            "frank_wolfe_penalized",
            "one_se_variance_retention",
            True,
            True,
            0.9,
        )
        self.assertEqual(report["minimum_loss_multiplier"], 64.0)
        self.assertEqual(report["best_multiplier"], 4.0)
        self.assertAlmostEqual(report["profile_variance_threshold"], 0.9)

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
