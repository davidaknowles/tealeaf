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
        self._assert_result(glm_solvers.fit_factorized(
            self.counts, self.compatibility, rank=2, max_iter=8, device="cpu", batch_cells=2
        ))

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

    def test_dense_admm(self):
        self._assert_result(glm_solvers.fit_admm(
            self.counts, self.compatibility, max_iter=8, inner_iter=5, device="cpu", max_dense_entries=100
        ))


if __name__ == "__main__":
    unittest.main()
