"""Small numerical checks for the Torch GLM backends."""

import unittest

import numpy as np
import scipy.sparse as sp

try:
    from tealeaf.sc import glm_solvers
    glm_solvers._torch()
except ImportError:  # pragma: no cover - exercised in non-GLM environments
    TORCH_AVAILABLE = False
else:
    TORCH_AVAILABLE = True


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
