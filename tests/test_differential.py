"""Tests for subisoform covariance and differential-splicing utilities."""

import tempfile
from pathlib import Path
import unittest

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import differential


class SpliceBlockTest(unittest.TestCase):
    def test_skipped_exon_becomes_two_path_block(self):
        lines = [
            'chr1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t1";',
            'chr1\tx\texon\t201\t300\t.\t+\t.\tgene_id "g"; transcript_id "t1";',
            'chr1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t2";',
            'chr1\tx\texon\t151\t170\t.\t+\t.\tgene_id "g"; transcript_id "t2";',
            'chr1\tx\texon\t201\t300\t.\t+\t.\tgene_id "g"; transcript_id "t2";',
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.gtf"
            path.write_text("\n".join(lines) + "\n")
            blocks = differential.build_splice_blocks(path)
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0].transcripts, ("t1", "t2"))
        self.assertEqual(blocks[0].n_paths, 2)
        self.assertNotEqual(*blocks[0].path_index)


class CovarianceTest(unittest.TestCase):
    def test_collapsing_recovers_identifiable_path_contrast(self):
        design = sp.csr_matrix(
            np.array([[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        )
        theta = np.array([0.2, 0.3, 0.5])
        path_index = np.array([0, 0, 1])
        result = differential.profiled_path_covariance(
            theta, path_index, (design,), (1000,)
        )
        self.assertTrue(result.identifiable)
        self.assertEqual(result.information_rank, 1)
        self.assertTrue(np.isfinite(result.covariance).all())
        self.assertGreater(result.covariance[0, 0], 0)

    def test_unidentifiable_path_contrast_is_reported(self):
        design = sp.csr_matrix(np.ones((1, 2)))
        result = differential.profiled_path_covariance(
            np.array([0.5, 0.5]),
            np.array([0, 1]),
            (design,),
            (1000,),
        )
        self.assertFalse(result.identifiable)
        self.assertTrue(np.isinf(result.covariance).all())

    def test_conditional_fisher_matches_parametric_variance(self):
        rng = np.random.default_rng(4)
        design = sp.eye(2, format="csr")
        baseline = np.array([0.3, 0.7])
        path_index = np.array([0, 1])
        total = 600
        analytic = differential.fit_path_perturbation(
            (np.array([180, 420]),),
            (design,),
            baseline,
            path_index,
        ).covariance.covariance[0, 0]
        self.assertAlmostEqual(
            analytic,
            1.0 / (2 * total * baseline[0] * baseline[1]),
            places=8,
        )
        estimates = []
        for _ in range(300):
            counts = rng.multinomial(total, baseline)
            fit = differential.fit_path_perturbation(
                (counts,), (design,), baseline, path_index
            )
            estimates.append(fit.path_logratios[0])
        empirical = np.var(estimates, ddof=1)
        self.assertAlmostEqual(empirical / analytic, 1.0, delta=0.25)

    def test_path_fit_moves_toward_observed_usage(self):
        design = sp.eye(2, format="csr")
        fit = differential.fit_path_perturbation(
            (np.array([80, 20]),),
            (design,),
            np.array([0.5, 0.5]),
            np.array([0, 1]),
        )
        self.assertTrue(fit.converged)
        np.testing.assert_allclose(fit.path_proportions, [0.8, 0.2], atol=1e-5)


class DifferentialTest(unittest.TestCase):
    def test_multivariate_gls_detects_condition_effect(self):
        rng = np.random.default_rng(7)
        condition = np.repeat([0, 1], 12)
        design = np.column_stack((np.ones(len(condition)), condition))
        values = rng.normal(0, 0.1, (len(condition), 2))
        values[condition == 1] += np.array([0.7, -0.4])
        covariance = np.repeat(
            (0.01 * np.eye(2))[None, :, :], len(condition), axis=0
        )
        result = differential.multivariate_gls_test(
            values, covariance, design, tested_columns=[1]
        )
        self.assertEqual(result["degrees_of_freedom"], 2)
        self.assertLess(result["p_value"], 1e-6)


if __name__ == "__main__":
    unittest.main()
