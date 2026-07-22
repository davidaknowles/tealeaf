"""Small numerical checks for the Torch GLM backends."""

import unittest
import gzip
from pathlib import Path
import tempfile

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import sc_utils

try:
    from tealeaf.sc import glm_solvers
    glm_solvers._torch()
except ImportError:  # pragma: no cover - exercised in non-GLM environments
    TORCH_AVAILABLE = False
else:
    TORCH_AVAILABLE = True


class ParameterizationTest(unittest.TestCase):
    def test_theta_design_is_length_rescaled_phi_design(self):
        membership = sp.csr_matrix(np.array([[1, 1], [0, 1]], dtype=float))
        inverse_lengths = np.array([0.5, 0.25])
        a_phi = sc_utils.weighted_ec_transcript_matrix(
            membership, inverse_lengths, parameterization="phi"
        )
        a_theta = sc_utils.weighted_ec_transcript_matrix(
            membership, inverse_lengths, parameterization="theta"
        )
        expected = a_phi @ sp.diags(1.0 / inverse_lengths)
        np.testing.assert_allclose(a_theta.toarray(), expected.toarray())

    def test_invalid_parameterization_is_rejected(self):
        with self.assertRaises(ValueError):
            sc_utils.weighted_ec_transcript_matrix(
                sp.eye(1), np.ones(1), parameterization="invalid"
            )

    def test_theta_nucnorm_uses_theta_design(self):
        counts = sp.csr_matrix(np.array([[8.0, 2.0], [2.0, 8.0]]))
        membership = sp.eye(2, format="csr")
        result = sc_utils.NNLS_nucnorm(
            counts,
            membership,
            np.array([0.5, 0.25]),
            regularization=0.0,
            max_iter=5,
            regularization_target="theta",
        )
        expected = counts.toarray() / np.asarray(counts.sum(axis=1))
        np.testing.assert_allclose(result, expected, atol=1e-7)

    def test_probability_sidecar_builds_column_normalized_design(self):
        membership = sp.csr_matrix(
            np.array([[1, 1], [1, 0], [0, 1]], dtype=float)
        )
        with tempfile.TemporaryDirectory() as directory:
            probability_file = Path(directory) / "probs.tsv.gz"
            cache_file = Path(directory) / "weights.npz"
            with gzip.open(probability_file, "wt") as handle:
                handle.write("cell_idx\teqid\tumi_rank\tprobs\n")
                handle.write("0\t0\t0\t8e-1,2e-1\n")
                handle.write("1\t0\t0\t6e-1,4e-1\n")
                handle.write("0\t1\t0\t1e0\n")
            design = sc_utils.averaged_ec_probability_matrix(
                probability_file, membership, cache_file
            )
            np.testing.assert_allclose(
                np.asarray(design.sum(axis=0)).ravel(), np.ones(2)
            )
            self.assertTrue(cache_file.is_file())
            cached = sc_utils.averaged_ec_probability_matrix(
                probability_file, membership, cache_file
            )
            np.testing.assert_allclose(cached.toarray(), design.toarray())

    def test_probability_sidecar_builds_separate_group_designs(self):
        membership = sp.csr_matrix(
            np.array([[1, 1], [1, 0], [0, 1]], dtype=float)
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            probability_file = directory / "probs.tsv.gz"
            caches = [directory / "poly.npz", directory / "hex.npz"]
            with gzip.open(probability_file, "wt") as handle:
                handle.write("cell_idx\teqid\tumi_rank\tprobs\n")
                handle.write("0\t0\t0\t9e-1,1e-1\n")
                handle.write("1\t0\t0\t2e-1,8e-1\n")
                handle.write("0\t1\t0\t1e0\n")
                handle.write("1\t2\t0\t1e0\n")
                handle.write("2\t0\t0\t5e-1,5e-1\n")
            designs = sc_utils.grouped_ec_probability_matrices(
                probability_file,
                membership,
                cell_groups=[0, 1, -1],
                n_groups=2,
                cache_files=caches,
            )
            self.assertEqual(len(designs), 2)
            for design in designs:
                np.testing.assert_allclose(
                    np.asarray(design.sum(axis=0)).ravel(), np.ones(2)
                )
            self.assertGreater(designs[0][0, 0], designs[1][0, 0])
            self.assertLess(designs[0][0, 1], designs[1][0, 1])
            self.assertTrue(all(path.is_file() for path in caches))


@unittest.skipUnless(TORCH_AVAILABLE, "Torch optional dependency is unavailable")
class GLMSolverTest(unittest.TestCase):
    def setUp(self):
        self.compatibility = sp.csr_matrix(
            np.array([[1.0, 0.2], [0.2, 1.0], [1.0, 1.0]], dtype=np.float32)
        )
        self.counts = sp.csr_matrix(
            np.array([[9, 1, 2], [1, 9, 2], [5, 5, 2], [0, 0, 0]], dtype=np.float32)
        )

    def _assert_result(self, result):
        output = glm_solvers.result_to_csr(result, 0, self.counts.shape[0]).toarray()
        self.assertTrue(np.isfinite(output).all())
        self.assertTrue((output >= -1e-7).all())
        np.testing.assert_allclose(output[:3].sum(axis=1), 1.0, atol=1e-4)
        np.testing.assert_allclose(output[3].sum(), 0.0, atol=1e-7)
        self.assertGreaterEqual(len(result.diagnostics["objective"]), 1)

    def test_sparse_glm_caches_normalized_blocks(self):
        data = glm_solvers.SparseGLM(
            self.counts, self.compatibility, device="cpu", batch_cells=2
        )
        first = list(data.blocks())
        second = list(data.blocks())
        self.assertIs(first[0][2], second[0][2])
        observed = np.vstack([block.to_dense().numpy() for _, _, block, _ in first])
        expected = self.counts.toarray()
        totals = expected.sum(axis=1, keepdims=True)
        expected = np.divide(
            expected, totals, out=np.zeros_like(expected), where=totals > 0
        )
        np.testing.assert_allclose(observed, expected, atol=1e-7)

    def test_sufficient_statistic_loss_matches_streamed_loss(self):
        torch = glm_solvers._torch()
        data = glm_solvers.SparseGLM(
            self.counts, self.compatibility, device="cpu", batch_cells=2
        )
        left = torch.tensor([[0.4, 0.1], [0.2, 0.3]])
        right = torch.tensor(
            [[0.2, 0.1], [0.3, 0.4], [0.1, 0.5], [0.0, 0.0]]
        )
        right_gram, ec_cross = data.factor_statistics(right)
        self.assertAlmostEqual(
            data.loss_for_factors(left, right),
            data.loss_from_statistics(left, right_gram, ec_cross),
            places=6,
        )

    def test_factorized_accepts_prepared_context_and_polishes(self):
        data = glm_solvers.SparseGLM(
            self.counts, self.compatibility, device="cpu", batch_cells=2
        )
        result = glm_solvers.fit_factorized(
            data,
            rank=2,
            max_iter=12,
            min_iter=12,
            polish_max_iter=3,
            device="cpu",
        )
        self.assertEqual(result.diagnostics["data_backend"], "cpu")
        self.assertIn("minibatch", result.diagnostics["phase"])
        self.assertEqual(result.diagnostics["phase"][-3:], ["polish"] * 3)
        self.assertGreater(result.diagnostics["cells_per_second"], 0)

    def test_factor_only_output(self):
        result = glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=2, max_iter=2,
            device="cpu", batch_cells=2,
        )
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "fit_"
            glm_solvers.write_chunked_result(
                result, prefix, ["a", "b", "c", "d"], ["t1", "t2"],
                write_chunks=False,
            )
            self.assertTrue(Path(f"{prefix}glm_factors.npz").is_file())
            self.assertFalse(list(Path(directory).glob("*glm_cells_*.npz")))
            self.assertGreater(result.diagnostics["active_cell_fraction"], 0)

    def test_factorized(self):
        result = glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=2, max_iter=8, device="cpu", batch_cells=2
        )
        self._assert_result(result)
        self.assertIsNone(result.diagnostics["factor_penalty"])
        left_norm = np.linalg.norm(result.left.numpy(), axis=0)
        right_norm = np.linalg.norm(result.right.numpy(), axis=0)
        active = (left_norm > 0) & (right_norm > 0)
        np.testing.assert_allclose(
            left_norm[active], right_norm[active], rtol=1e-5, atol=1e-7
        )

    def test_factorized_expands_a_rank_warm_start(self):
        initial = glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=1, max_iter=2,
            min_iter=2, device="cpu", batch_cells=2,
        )
        expanded = glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=3, max_iter=2,
            min_iter=2, device="cpu", batch_cells=2,
            initial_factors=(initial.left, initial.right),
        )
        self.assertTrue(expanded.diagnostics["warm_started"])
        self.assertEqual(expanded.diagnostics["warm_start_rank"], 1)
        self.assertEqual(expanded.left.shape[1], 3)

    def test_factorized_uses_minimum_iterations_and_patience(self):
        result = glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=2, max_iter=20,
            min_iter=4, patience=2, tol=1.0, device="cpu", batch_cells=2,
        )
        self.assertGreaterEqual(result.diagnostics["iterations"], 4)
        self.assertTrue(result.diagnostics["converged"])
        self.assertEqual(
            result.diagnostics["convergence_reason"], "objective_patience"
        )

    def test_factorized_admm(self):
        self._assert_result(glm_solvers.fit_factorized_admm(
            self.counts, self.compatibility, rank=2, max_iter=8, device="cpu", batch_cells=2
        ))

    def test_factorized_admm_adapts_rho_and_rescales_duals(self):
        result = glm_solvers.fit_factorized_admm(
            self.counts,
            self.compatibility,
            rank=2,
            rho=1e-4,
            max_iter=4,
            min_iter=4,
            rho_update_interval=1,
            rho_balance=1.01,
            rho_scale=2.0,
            device="cpu",
            batch_cells=2,
        )
        self.assertTrue(result.diagnostics["adaptive_rho"])
        self.assertGreater(len(result.diagnostics["rho_updates"]), 0)
        self.assertNotEqual(result.diagnostics["final_rho"], 1e-4)
        self.assertTrue(np.isfinite(result.diagnostics["final_rho"]))

    def test_factorized_admm_accepts_primal_warm_start(self):
        initial = glm_solvers.fit_factorized_admm(
            self.counts, self.compatibility, rank=2, max_iter=2,
            min_iter=2, device="cpu", batch_cells=2,
        )
        continued = glm_solvers.fit_factorized_admm(
            self.counts, self.compatibility, rank=2, max_iter=2,
            min_iter=2, device="cpu", batch_cells=2,
            initial_factors=(initial.left, initial.right),
        )
        self.assertTrue(continued.diagnostics["warm_started"])
        self.assertEqual(continued.diagnostics["warm_start_rank"], 2)

    def test_factor_loss_ignores_empty_response_cells(self):
        torch = glm_solvers._torch()
        counts = self.counts.copy().tolil()
        counts[0, :] = 0
        data = glm_solvers.SparseGLM(
            counts.tocsr(), self.compatibility, device="cpu", batch_cells=2
        )
        left = torch.tensor([[0.4], [0.2]])
        right = torch.tensor([[100.0], [0.3], [0.2], [0.1]])
        expected_right = right.clone()
        expected_right[0] = 0
        self.assertAlmostEqual(
            data.loss_for_factors(left, right),
            data.loss_for_factors(left, expected_right),
            places=6,
        )

    def test_profile_diagnostics_reject_library_size_only_variation(self):
        torch = glm_solvers._torch()
        collapsed = glm_solvers.GLMResult(
            "test",
            torch.tensor([[1.0], [2.0]]),
            torch.tensor([[1.0], [2.0], [3.0]]),
            None,
            {},
        )
        collapsed_stats = glm_solvers.factor_profile_diagnostics(collapsed)
        self.assertEqual(collapsed_stats["normalized_profile_active_fraction"], 1.0)
        self.assertLess(
            collapsed_stats["normalized_profile_relative_variance"], 1e-6
        )
        variable = glm_solvers.GLMResult(
            "test",
            torch.eye(2),
            torch.tensor([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
            None,
            {},
        )
        variable_stats = glm_solvers.factor_profile_diagnostics(variable)
        self.assertGreater(
            variable_stats["normalized_profile_relative_variance"], 1e-3
        )

    def test_frank_wolfe(self):
        result = glm_solvers.fit_frank_wolfe(
            self.counts, self.compatibility, rank=3, max_iter=3, device="cpu", batch_cells=2
        )
        self._assert_result(result)
        nuclear_norm = np.linalg.svd(
            result.left.detach().cpu().numpy() @ result.right.detach().cpu().numpy().T,
            compute_uv=False,
        ).sum()
        self.assertLessEqual(nuclear_norm, result.diagnostics["tau"] + 1e-4)

    def test_frank_wolfe_uses_duality_gap_patience(self):
        result = glm_solvers.fit_frank_wolfe(
            self.counts, self.compatibility, rank=2, max_atoms=10,
            max_iter=10, min_iter=4, patience=2, tol=1e6,
            device="cpu", batch_cells=2,
        )
        self.assertGreaterEqual(result.diagnostics["iterations"], 4)
        self.assertTrue(result.diagnostics["converged"])
        self.assertEqual(len(result.diagnostics["duality_gap"]), 4)
        self.assertEqual(
            result.diagnostics["convergence_reason"], "duality_gap_patience"
        )

    def test_penalized_frank_wolfe_uses_signed_nuclear_oracle(self):
        result = glm_solvers.fit_frank_wolfe_penalized(
            self.counts,
            self.compatibility,
            rank=4,
            max_atoms=4,
            max_iter=4,
            power_iter=8,
            nonnegative_penalty=1.0,
            device="cpu",
            batch_cells=2,
        )
        self._assert_result(result)
        fitted = (
            result.left.detach().cpu().numpy()
            @ result.right.detach().cpu().numpy().T
        )
        self.assertLessEqual(
            np.linalg.svd(fitted, compute_uv=False).sum(),
            result.diagnostics["tau"] + 1e-4,
        )
        self.assertFalse(result.diagnostics["gap_is_certificate"])
        self.assertEqual(len(result.diagnostics["candidate_gap"]), 4)
        self.assertGreaterEqual(
            result.diagnostics["negative_l1_fraction"], 0.0
        )
        self.assertLessEqual(
            result.diagnostics["negative_l1_fraction"], 1.0
        )

    def test_penalized_frank_wolfe_stops_on_objective_patience(self):
        result = glm_solvers.fit_frank_wolfe_penalized(
            self.counts, self.compatibility, rank=2, max_atoms=10,
            max_iter=10, min_iter=4, patience=2, tol=1e6,
            device="cpu", batch_cells=2,
        )
        self.assertTrue(result.diagnostics["converged"])
        self.assertEqual(result.diagnostics["iterations"], 4)
        self.assertEqual(
            result.diagnostics["convergence_reason"], "objective_patience"
        )

    def test_penalized_frank_wolfe_continues_prior_atoms(self):
        initial = glm_solvers.fit_frank_wolfe_penalized(
            self.counts, self.compatibility, tau=0.5, max_atoms=6,
            max_iter=2, min_iter=2, device="cpu", batch_cells=2,
        )
        continued = glm_solvers.fit_frank_wolfe_penalized(
            self.counts, self.compatibility, tau=1.0, max_atoms=6,
            max_iter=2, min_iter=2, device="cpu", batch_cells=2,
            initial_factors=(initial.left, initial.right),
        )
        self.assertTrue(continued.diagnostics["warm_started"])
        self.assertGreater(continued.diagnostics["warm_start_rank"], 0)
        self.assertGreater(
            continued.left.shape[1], continued.diagnostics["warm_start_rank"]
        )

    def test_penalized_frank_wolfe_operator_matches_dense_gradient(self):
        torch = glm_solvers._torch()
        data = glm_solvers.SparseGLM(
            self.counts, self.compatibility, device="cpu", batch_cells=2
        )
        left = torch.tensor([[0.4, -0.2], [0.1, 0.3]])
        right = torch.tensor(
            [[0.2, -0.1], [-0.3, 0.5], [0.4, 0.2], [0.0, 0.0]]
        )
        penalty = 0.7
        normalized = self.counts.toarray()
        totals = normalized.sum(axis=1, keepdims=True)
        normalized = np.divide(
            normalized, totals, out=np.zeros_like(normalized), where=totals > 0
        )
        compatibility = self.compatibility.toarray()
        phi = left.numpy() @ right.numpy().T
        q_dense = (
            compatibility.T @ normalized.T
            - compatibility.T @ compatibility @ phi
            - penalty * np.minimum(phi, 0.0)
        )
        cell_vectors = torch.tensor(
            [[0.1, 0.2], [-0.2, 0.4], [0.3, -0.1], [0.5, 0.2]]
        )
        transcript_vectors = torch.tensor([[0.3, -0.1], [0.2, 0.4]])
        np.testing.assert_allclose(
            glm_solvers._penalized_fw_q_times(
                data, left, right, cell_vectors, penalty
            ).numpy(),
            q_dense @ cell_vectors.numpy(),
            atol=1e-6,
        )
        np.testing.assert_allclose(
            glm_solvers._penalized_fw_qt_times(
                data, left, right, transcript_vectors, penalty
            ).numpy(),
            q_dense.T @ transcript_vectors.numpy(),
            atol=1e-6,
        )

    def test_invalid_stopping_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            glm_solvers.fit_factorized(
                self.counts, self.compatibility, max_iter=0, device="cpu"
            )

    def test_factorized_rejects_rank_reduction_warm_start(self):
        initial = glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=2, max_iter=1,
            min_iter=1, device="cpu",
        )
        with self.assertRaisesRegex(ValueError, "exceeds requested rank"):
            glm_solvers.fit_factorized(
                self.counts, self.compatibility, rank=1, max_iter=1,
                min_iter=1, device="cpu",
                initial_factors=(initial.left, initial.right),
            )

    def test_dense_admm(self):
        self._assert_result(glm_solvers.fit_admm(
            self.counts, self.compatibility, max_iter=8, inner_iter=5, device="cpu", max_dense_entries=100
        ))


if __name__ == "__main__":
    unittest.main()
