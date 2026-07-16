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

    def test_dense_admm(self):
        self._assert_result(glm_solvers.fit_admm(
            self.counts, self.compatibility, max_iter=8, inner_iter=5, device="cpu", max_dense_entries=100
        ))


if __name__ == "__main__":
    unittest.main()
