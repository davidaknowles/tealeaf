"""Checks for scalable GLM count-fold cross-validation."""

import unittest
from unittest import mock
from pathlib import Path
import tempfile

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import glm_cv

try:
    from tealeaf.sc import glm_solvers
    glm_solvers._torch()
except ImportError:  # pragma: no cover
    TORCH_AVAILABLE = False
else:
    TORCH_AVAILABLE = True


class PairedPrimerPreparationTest(unittest.TestCase):
    def test_sparse_storage_bytes_counts_csr_arrays(self):
        matrix = sp.csr_matrix([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32)
        expected = (
            matrix.data.nbytes
            + matrix.indices.nbytes
            + matrix.indptr.nbytes
        )
        self.assertEqual(glm_cv.sparse_storage_bytes(matrix), expected)

    def test_paired_response_and_design_use_equal_primer_weight(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            membership = sp.eye(2, format="csr")
            counts = sp.csr_matrix(np.array(
                [[8, 2], [2, 8], [4, 0], [0, 3]], dtype=np.int64
            ))
            sp.save_npz(directory / "gene_eqclass.npz", membership)
            sp.save_npz(directory / "geqc_counts.npz", counts)
            (directory / "quants_mat_cols.txt").write_text("tx1\ntx2\n")
            (directory / "quants_mat_rows.txt").write_text(
                "poly1\nhex1\npoly2\nhex2\n"
            )
            fasta = directory / "transcripts.fa"
            fasta.write_text(">tx1\n" + "A" * 400 + "\n>tx2\n" + "C" * 400 + "\n")
            pairs = directory / "pairs.tsv"
            pairs.write_text(
                "cell_id\tpolydt_barcode\tranhex_barcode\n"
                "cell1\tpoly1\thex1\n"
                "cell1_duplicate\tpoly1\thex1\n"
                "cell2\tpoly2\thex2\n"
            )
            prepared = glm_cv.prepare_paired_primer_glm_data(
                directory,
                fasta,
                pairs,
                ec_design="binary",
                regularization_target="phi",
                min_eq=1,
                min_half_umis=5,
            )
            self.assertEqual(prepared.counts.shape, (1, 4))
            np.testing.assert_allclose(
                prepared.counts.toarray(), [[0.4, 0.1, 0.1, 0.4]]
            )
            np.testing.assert_allclose(
                prepared.compatibility.toarray(),
                [[0.5, 0], [0, 0.5], [0.5, 0], [0, 0.5]],
            )
            np.testing.assert_array_equal(prepared.barcodes, ["cell1"])
            self.assertEqual(prepared.metadata["retained_pair_count"], 1)
            np.testing.assert_array_equal(
                prepared.cv_raw_counts.toarray(),
                [[8, 2, 2, 8]],
            )

    def test_paired_count_folds_split_raw_molecules_before_normalizing(self):
        raw = sp.csr_matrix(
            np.array([[8, 2, 2, 8], [4, 6, 7, 3]], dtype=np.int64)
        )
        pairs = glm_cv.paired_primer_count_fold_pairs(
            raw, n_folds=2, seed=3
        )
        self.assertEqual(len(pairs), 2)
        for training, validation in pairs:
            np.testing.assert_allclose(
                np.asarray(training[:, :2].sum(axis=1)).ravel(), 0.5
            )
            np.testing.assert_allclose(
                np.asarray(training[:, 2:].sum(axis=1)).ravel(), 0.5
            )
            np.testing.assert_allclose(
                np.asarray(validation[:, :2].sum(axis=1)).ravel(), 0.5
            )
            np.testing.assert_allclose(
                np.asarray(validation[:, 2:].sum(axis=1)).ravel(), 0.5
            )

    def test_positional_design_loads_separate_primer_caches(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            membership = sp.csr_matrix(np.ones((1, 2)))
            sp.save_npz(directory / "gene_eqclass.npz", membership)
            sp.save_npz(
                directory / "geqc_counts.npz",
                sp.csr_matrix(np.array([[10], [10]])),
            )
            (directory / "quants_mat_cols.txt").write_text("tx1\ntx2\n")
            (directory / "quants_mat_rows.txt").write_text("poly\nhex\n")
            (directory / "transcripts.fa").write_text(
                ">tx1\n" + "A" * 400 + "\n>tx2\n" + "C" * 400 + "\n"
            )
            (directory / "pairs.tsv").write_text(
                "cell_id\tpolydt_barcode\tranhex_barcode\ncell\tpoly\thex\n"
            )
            sp.save_npz(
                directory / "gene_eqclass_posbias_polydt.npz",
                sp.csr_matrix([[0.8, 0.2]]),
            )
            sp.save_npz(
                directory / "gene_eqclass_posbias_ranhex.npz",
                sp.csr_matrix([[0.1, 0.9]]),
            )
            np.save(
                directory / "salmon_effective_lengths_polydt.npy",
                np.array([100.0, 200.0]),
            )
            np.save(
                directory / "salmon_effective_lengths_ranhex.npy",
                np.array([150.0, 250.0]),
            )
            prepared = glm_cv.prepare_paired_primer_glm_data(
                directory,
                directory / "transcripts.fa",
                directory / "pairs.tsv",
                ec_design="positional",
                regularization_target="phi",
                min_eq=1,
                min_half_umis=1,
            )
            np.testing.assert_allclose(
                prepared.compatibility.toarray(),
                [[0.5, 0.5], [0.5, 0.5]],
            )
            self.assertEqual(prepared.metadata["ec_design"], "positional")
            theta = glm_cv.prepare_paired_primer_glm_data(
                directory,
                directory / "transcripts.fa",
                directory / "pairs.tsv",
                ec_design="positional",
                regularization_target="theta",
                primer_sampling_model="oligodt_tpm",
                min_eq=1,
                min_half_umis=1,
            )
            np.testing.assert_allclose(
                theta.compatibility.toarray(),
                [[0.5, 0.5], [0.375, 0.625]],
            )
            self.assertEqual(
                theta.metadata["primer_sampling_model"], "oligodt_tpm"
            )
            self.assertEqual(
                theta.metadata["primer_sampling_scales"], [1.0, 200.0]
            )


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

    def test_cell_sampling_uses_supplied_raw_count_threshold(self):
        raw_totals = np.array([10, 500, 499, 900])
        selected = glm_cv.sample_cells_by_count(
            self.counts,
            0,
            min_count=500,
            totals=raw_totals,
        )
        np.testing.assert_array_equal(selected, [1, 3])

    def test_cell_sampling_validates_supplied_totals(self):
        with self.assertRaisesRegex(ValueError, "one value per"):
            glm_cv.sample_cells_by_count(
                self.counts, 0, min_count=1, totals=[1, 2]
            )

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

    def test_adaptive_grid_replays_warm_start_path_after_expansion(self):
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
            "multipliers": [1.0, 2.0, 8.0],
            "fold_scales": [1.0, 1.0],
            "mean_validation_loss": {1.0: 0.2, 2.0: 0.1, 8.0: 0.3},
            "validation_standard_error": {1.0: 0.01, 2.0: 0.01, 8.0: 0.01},
            "candidate_converged": {1.0: True, 2.0: True, 8.0: True},
            "candidate_nondegenerate": {1.0: True, 2.0: True, 8.0: True},
            "mean_profile_relative_variance": {1.0: 0.1, 2.0: 0.2, 8.0: 0.3},
            "minimum_profile_active_fraction": {1.0: 1.0, 2.0: 1.0, 8.0: 1.0},
            "best_multiplier": 2.0,
            "best_on_boundary": False,
            "fold_results": [
                {"multiplier": 1.0}, {"multiplier": 2.0},
                {"multiplier": 8.0},
            ],
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
        self.assertEqual(mocked.call_args_list[1].args[3], [1.0, 2.0, 8.0])
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
                "rank": 4,
                "max_atoms": 4,
                "max_iter": 2,
                "min_iter": 2,
                "power_iter": 3,
            },
        )
        self.assertIn(report["best_multiplier"], [0.5, 2.0])
        self.assertIsInstance(report["best_on_boundary"], bool)
        self.assertEqual(len(report["fold_results"]), 4)
        self.assertEqual(report["regularization_path"], [0.5, 2.0])
        for fold in range(2):
            rows = [row for row in report["fold_results"] if row["fold"] == fold]
            self.assertEqual([row["multiplier"] for row in rows], [0.5, 2.0])
            self.assertFalse(rows[0]["warm_started"])
            self.assertTrue(rows[1]["warm_started"])

    def test_admm_path_runs_from_strong_to_weak_regularization(self):
        report = glm_cv.cross_validate_glm(
            self.counts,
            self.compatibility,
            "admm_factorized",
            [0.0, 1e-4],
            n_folds=2,
            device="cpu",
            batch_cells=2,
            power_iter=2,
            fit_kwargs={
                "rank": 2,
                "max_iter": 2,
                "min_iter": 2,
            },
        )
        self.assertEqual(report["regularization_path"], [1e-4, 0.0])
        for fold in range(2):
            rows = [row for row in report["fold_results"] if row["fold"] == fold]
            self.assertEqual([row["multiplier"] for row in rows], [1e-4, 0.0])
            self.assertFalse(rows[0]["warm_started"])
            self.assertTrue(rows[1]["warm_started"])

    def test_factorized_rank_cv_runs_from_few_to_many_factors(self):
        progress = []
        report = glm_cv.cross_validate_factorized_rank(
            self.counts,
            self.compatibility,
            [1, 2],
            n_folds=2,
            device="cpu",
            batch_cells=2,
            fit_kwargs={"max_iter": 2, "min_iter": 2},
            selection_rule="one_standard_error",
            progress_callback=progress.append,
        )
        self.assertEqual(report["rank_path"], [1, 2])
        self.assertIn(report["best_rank"], [1, 2])
        self.assertIsNone(report["regularization"])
        for fold in range(2):
            rows = [row for row in report["fold_results"] if row["fold"] == fold]
            self.assertEqual([row["rank"] for row in rows], [1, 2])
            self.assertFalse(rows[0]["warm_started"])
            self.assertTrue(rows[1]["warm_started"])
            self.assertEqual(rows[1]["warm_start_rank"], 1)
        self.assertEqual(
            [(row["fold"], row["rank"]) for row in progress],
            [(0, 1), (0, 2), (1, 1), (1, 2)],
        )

    def test_adaptive_rank_cv_continues_without_replaying_completed_ranks(self):
        initial = {"best_rank": 2}
        extended = {"best_rank": 2}
        state = {"continuation": True}
        with mock.patch.object(
            glm_cv,
            "cross_validate_factorized_rank",
            side_effect=[(initial, state), (extended, state)],
        ) as mocked:
            report = glm_cv.cross_validate_factorized_rank_adaptive(
                self.counts,
                self.compatibility,
                [1, 2],
                max_rank=8,
                max_grid_expansions=2,
            )
        self.assertEqual(mocked.call_count, 2)
        self.assertEqual(mocked.call_args_list[1].args[2], [1, 2, 4])
        self.assertIs(mocked.call_args_list[1].kwargs["_state"], state)
        self.assertEqual(report["grid_expansions"], 1)


if __name__ == "__main__":
    unittest.main()
