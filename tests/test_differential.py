"""Tests for subisoform covariance and differential-splicing utilities."""

import tempfile
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
import scipy.special
import scipy.sparse as sp

from extra_scripts import (
    merge_celltype_compositional_splicing,
    run_celltype_compositional_splicing,
    run_differential_splicing,
)
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

    def test_shared_cell_fit_matches_homogeneous_pseudobulk(self):
        design = sp.eye(2, format="csr")
        cell_counts = sp.csr_matrix(
            np.array([[30, 70], [20, 80], [40, 60]])
        )
        baseline = np.repeat([[0.5, 0.5]], 3, axis=0)
        path_index = np.array([0, 1])
        shared = differential.fit_shared_path_perturbation(
            (cell_counts,),
            (design,),
            baseline,
            path_index,
            weights=np.ones(3),
        )
        aggregated = differential.fit_path_perturbation(
            (np.asarray(cell_counts.sum(axis=0)).ravel(),),
            (design,),
            np.array([0.5, 0.5]),
            path_index,
        )
        self.assertTrue(shared.converged)
        np.testing.assert_allclose(
            shared.path_proportions,
            aggregated.path_proportions,
            atol=1e-6,
        )
        np.testing.assert_allclose(
            shared.covariance.covariance,
            aggregated.covariance.covariance,
            rtol=1e-5,
        )

    def test_shared_cell_fit_recovers_shift_with_distinct_baselines(self):
        design = sp.eye(2, format="csr")
        baselines = np.array([[0.8, 0.2], [0.2, 0.8]])
        path_index = np.array([0, 1])
        basis = differential.helmert_basis(2)
        shifted, _ = differential._perturbed_theta_matrix(
            baselines,
            path_index,
            basis,
            np.array([0.7]),
        )
        counts = sp.csr_matrix(np.rint(100_000 * shifted).astype(int))
        fit = differential.fit_shared_path_perturbation(
            (counts,),
            (design,),
            baselines,
            path_index,
            weights=np.ones(2),
        )
        self.assertTrue(fit.converged)
        np.testing.assert_allclose(fit.delta, [0.7], atol=1e-4)


class ClusteredCompositionalTest(unittest.TestCase):
    @staticmethod
    def simulated_counts(seed=17, subjects=36, total=80, effect=0.9):
        rng = np.random.default_rng(seed)
        clusters = np.repeat(np.arange(subjects), 2)
        condition = np.repeat(np.arange(subjects) >= subjects // 2, 2)
        cell_type = np.tile([0, 1], subjects)
        null_design = np.c_[1 - cell_type, cell_type]
        alternative_design = np.c_[null_design, condition]
        random_intercepts = rng.normal(0.0, 0.55, size=subjects)
        logits = (
            -0.4 * null_design[:, 0]
            + 0.3 * null_design[:, 1]
            + effect * condition
            + random_intercepts[clusters]
        )
        probabilities = scipy.special.expit(logits)
        successes = rng.binomial(total, probabilities)
        counts = np.c_[successes, total - successes].astype(float)
        return counts, null_design, alternative_design, clusters

    def test_multinomial_gee_recovers_clustered_effect(self):
        counts, _, design, clusters = self.simulated_counts()
        result = differential.multinomial_gee_test(
            counts,
            design,
            tested_columns=[2],
            clusters=clusters,
        )
        self.assertTrue(result["converged"])
        self.assertGreater(result["working_correlation"], 0)
        self.assertAlmostEqual(result["coefficients"][2, 0], 0.9, delta=0.35)
        self.assertLess(result["p_value"], 0.05)
        warm = differential.multinomial_gee_test(
            counts,
            design,
            tested_columns=[2],
            clusters=clusters,
            initial=result["coefficients"],
        )
        self.assertAlmostEqual(warm["statistic"], result["statistic"], places=5)

    def test_paired_logratio_recovers_subject_level_shift(self):
        rng = np.random.default_rng(41)
        rows = []
        labels = []
        clusters = []
        for subject in range(12):
            baseline = rng.normal(0.0, 0.3)
            for label, effect in ((0, 0.0), (1, 0.8)):
                probability = scipy.special.expit(baseline + effect)
                first = rng.binomial(100, probability)
                rows.append([first, 100 - first])
                labels.append(label)
                clusters.append(subject)
        result = differential.paired_logratio_test(rows, labels, clusters)
        self.assertTrue(result["converged"])
        self.assertEqual(result["n_subjects"], 12)
        self.assertLess(result["p_value"], 0.01)

    def test_multinomial_glmm_recovers_clustered_effect(self):
        counts, null, alternative, clusters = self.simulated_counts(
            subjects=28,
            total=60,
        )
        result = differential.multinomial_glmm_test(
            counts,
            null,
            alternative,
            clusters,
            max_iter=100,
        )
        self.assertTrue(result["null_converged"])
        self.assertTrue(result["alternative_converged"])
        self.assertGreater(result["null_random_effect_sd"], 0.1)
        self.assertLess(result["p_value"], 0.05)

    def test_integerized_counts_preserve_rounded_totals(self):
        counts = np.array([[1.2, 2.2, 3.2], [0.1, 0.1, 0.1]])
        rounded = differential.integerize_compositional_counts(counts)
        np.testing.assert_array_equal(rounded.sum(axis=1), [7, 1])
        self.assertTrue(np.issubdtype(rounded.dtype, np.integer))
        self.assertTrue((rounded >= 0).all())

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

    def test_clustered_multivariate_wald_detects_paired_effect(self):
        rng = np.random.default_rng(31)
        clusters = np.repeat(np.arange(30), 2)
        tested = np.tile([0, 1], 30)
        values = rng.normal(0, 0.15, (len(clusters), 2))
        values += rng.normal(0, 0.3, (30, 2))[clusters]
        values[tested == 1] += np.array([0.6, -0.4])
        covariances = np.repeat(
            (0.02 * np.eye(2))[None, :, :], len(clusters), axis=0
        )
        result = differential.clustered_multivariate_wald_test(
            values,
            covariances,
            np.column_stack((np.ones(len(clusters)), tested)),
            tested_columns=[1],
            clusters=clusters,
        )
        self.assertEqual(result["degrees_of_freedom"], 2)
        self.assertLess(result["p_value"], 1e-6)

    def test_clustered_multivariate_wald_leverage_adjustments_run(self):
        rng = np.random.default_rng(33)
        clusters = np.repeat(np.arange(12), 2)
        tested = np.tile([0, 1], 12)
        values = rng.normal(size=(len(clusters), 2))
        covariances = np.repeat(np.eye(2)[None, :, :], len(values), axis=0)
        design = np.column_stack((np.ones(len(values)), tested))
        for adjustment in ("cr2", "cr3"):
            result = differential.clustered_multivariate_wald_test(
                values,
                covariances,
                design,
                tested_columns=[1],
                clusters=clusters,
                cluster_adjustment=adjustment,
            )
            self.assertEqual(result["cluster_adjustment"], adjustment)
            self.assertTrue(np.isfinite(result["p_value"]))

    def test_clustered_multivariate_wald_requires_residual_cluster_df(self):
        values = np.arange(8, dtype=float).reshape(4, 2)
        covariances = np.repeat(np.eye(2)[None, :, :], 4, axis=0)
        design = np.column_stack((np.ones(4), [0, 0, 1, 1]))
        with self.assertRaisesRegex(ValueError, "cluster count"):
            differential.clustered_multivariate_wald_test(
                values,
                covariances,
                design,
                tested_columns=[1],
                clusters=[0, 0, 1, 1],
            )

    def test_restricted_wild_cluster_null_preserves_cluster_residuals(self):
        clusters = np.repeat(np.arange(6), 2)
        nuisance = np.ones((len(clusters), 1))
        values = np.column_stack((clusters, -clusters)).astype(float)
        covariances = np.repeat(np.eye(2)[None, :, :], len(values), axis=0)
        simulated = differential.restricted_wild_cluster_null(
            values,
            covariances,
            nuisance,
            clusters,
            np.random.default_rng(2),
        )
        fitted = np.repeat(values.mean(axis=0)[None, :], len(values), axis=0)
        np.testing.assert_allclose(
            np.abs(simulated - fitted), np.abs(values - fitted)
        )

    def test_effective_multinomial_size_recovers_known_size(self):
        proportions = np.array([0.2, 0.3, 0.5])
        basis = differential.helmert_basis(3)
        covariance = (
            basis.T @ np.diag(1.0 / proportions) @ basis / 250.0
        )
        size = differential.effective_multinomial_size(
            proportions, covariance
        )
        self.assertAlmostEqual(size, 250.0)

    def test_composition_projection_preserves_selected_logratio(self):
        proportions = np.array([0.2, 0.3, 0.5])
        covariance = np.array([[0.4, 0.1], [0.1, 0.6]])
        projected, projected_covariance = (
            differential.project_composition(
                proportions, covariance, [0, 2]
            )
        )
        np.testing.assert_allclose(projected, [2 / 7, 5 / 7])
        old_basis = differential.helmert_basis(3)
        contrast = np.array([1.0, 0.0, -1.0]) / np.sqrt(2)
        expected = (
            contrast @ old_basis @ covariance
            @ old_basis.T @ contrast
        )
        self.assertAlmostEqual(projected_covariance[0, 0], expected)

    def test_permutation_rank_treats_optimizer_noise_as_ties(self):
        p_value = differential.permutation_rank_p_value(
            2.0 + 1e-8,
            [0.5, 2.0, 2.0 - 1e-8],
        )
        self.assertEqual(p_value, 0.75)

    def test_empirical_null_calibration_excludes_current_block(self):
        table = pd.DataFrame({
            "block_id": ["b1", "b2"],
            "degrees_of_freedom": [1, 1],
            "asymptotic_p_value": [0.01, 0.03],
        })
        null = pd.DataFrame({
            "block_id": ["b1", "b1", "b2", "b2"],
            "degrees_of_freedom": [1, 1, 1, 1],
            "p_value": [0.01, 0.20, 0.02, 0.30],
        })
        calibrated, null_calibrated = (
            merge_celltype_compositional_splicing
            .empirical_null_calibration(table, null)
        )
        np.testing.assert_allclose(
            calibrated["pooled_permutation_p_value"],
            [1 / 3, 2 / 3],
        )
        self.assertEqual(len(null_calibrated), 4)

    def test_pairwise_design_requires_shared_subjects(self):
        records = [
            {
                "cluster": cell_type,
                "condition": "control",
                "mouse": mouse,
            }
            for mouse in ("m1", "m2", "m3")
            for cell_type in ("A", "B")
        ]
        records.append({
            "cluster": "C",
            "condition": "control",
            "mouse": "m4",
        })
        designs = (
            run_celltype_compositional_splicing
            .pairwise_celltype_designs(records, min_mice=3)
        )
        self.assertEqual(len(designs), 1)
        contrast, prepared = designs[0]
        self.assertEqual(contrast, "A_vs_B")
        self.assertEqual(prepared[0], records[:6])
        self.assertEqual(prepared[4].shape, (6, 3))
        self.assertEqual(prepared[5].shape, (6, 4))

    def test_dirichlet_multinomial_detects_composition_effect(self):
        counts = np.array([
            [85, 15],
            [81, 19],
            [88, 12],
            [15, 85],
            [21, 79],
            [12, 88],
        ])
        alternative = np.c_[
            np.ones(6),
            np.r_[np.zeros(3), np.ones(3)],
        ]
        result = differential.dirichlet_multinomial_test(
            counts,
            np.ones((6, 1)),
            alternative,
        )
        self.assertTrue(result["null_converged"])
        self.assertTrue(result["alternative_converged"])
        self.assertLess(result["p_value"], 0.05)
        reused = differential.dirichlet_multinomial_test(
            counts,
            np.ones((6, 1)),
            alternative,
            fitted_null=result,
        )
        self.assertAlmostEqual(reused["statistic"], result["statistic"])
        free = differential.dirichlet_multinomial_test(
            counts,
            np.ones((6, 1)),
            alternative,
            fix_null_concentration=False,
        )
        self.assertEqual(
            free["degrees_of_freedom"], result["degrees_of_freedom"]
        )
        self.assertGreaterEqual(free["statistic"], result["statistic"])
        free_reused = differential.dirichlet_multinomial_test(
            counts,
            np.ones((6, 1)),
            alternative,
            fix_null_concentration=False,
            fitted_null=free,
        )
        self.assertAlmostEqual(
            free_reused["statistic"], free["statistic"]
        )

    def test_clustered_gls_detects_subject_level_effect(self):
        rng = np.random.default_rng(21)
        subjects = np.repeat(np.arange(20), 2)
        cell_type = np.tile([0, 1], 20)
        condition = np.repeat(np.r_[np.zeros(10), np.ones(10)], 2)
        subject_effect = np.repeat(rng.normal(0, 0.3, 20), 2)
        values = (
            0.4 * cell_type
            + 0.8 * condition
            + subject_effect
            + rng.normal(0, 0.05, 40)
        )[:, None]
        design = np.column_stack((
            cell_type == 0,
            cell_type == 1,
            condition,
        ))
        covariance = np.repeat(
            np.array([[[0.01]]]), len(values), axis=0
        )
        result = differential.clustered_multivariate_gls_test(
            values,
            covariance,
            design,
            tested_columns=[2],
            clusters=subjects,
        )
        self.assertGreater(result["cluster_variance"], 0.01)
        self.assertLess(result["p_value"], 0.01)
        fixed = differential.clustered_multivariate_gls_test(
            values,
            covariance,
            design,
            tested_columns=[2],
            clusters=subjects,
            variance_components=(
                result["cluster_variance"],
                result["residual_variance"],
            ),
        )
        self.assertAlmostEqual(fixed["statistic"], result["statistic"])

    def test_gaussian_mixed_model_lrt_and_score_bf_detect_effect(self):
        rng = np.random.default_rng(41)
        subjects = 24
        clusters = np.repeat(np.arange(subjects), 2)
        treatment = np.tile([0.0, 1.0], subjects)
        subject_effect = np.repeat(rng.normal(0, 0.25, subjects), 2)
        values = (
            0.8 * treatment
            + subject_effect
            + rng.normal(0, 0.08, len(treatment))
        )[:, None]
        covariances = np.repeat(
            np.array([[[0.01]]]), len(values), axis=0
        )
        result = differential.gaussian_mixed_model_lrt(
            values,
            covariances,
            np.column_stack((np.ones(len(values)), treatment)),
            tested_columns=[1],
            clusters=clusters,
            effect_prior_scales=(0.25, 0.5, 1.0),
        )
        self.assertLess(result["p_value"], 0.01)
        self.assertGreater(result["mixture_log_bayes_factor"], np.log(10.0))
        self.assertEqual(result["degrees_of_freedom"], 1)
        self.assertEqual(result["efficient_score"].shape, (1,))
        self.assertEqual(result["efficient_information"].shape, (1, 1))
        self.assertGreater(result["alternative_cluster_variance"], 0.01)

        def dense_profile_objective(model_design, cluster_variance, residual_variance):
            precision = np.zeros((model_design.shape[1], model_design.shape[1]))
            score = np.zeros(model_design.shape[1])
            quadratic = 0.0
            log_determinant = 0.0
            for cluster in np.unique(clusters):
                rows = np.flatnonzero(clusters == cluster)
                variance = np.diag(
                    covariances[rows, 0, 0] + residual_variance
                ) + cluster_variance * np.ones((len(rows), len(rows)))
                weight = np.linalg.inv(variance)
                local_values = values[rows, 0]
                precision += model_design[rows].T @ weight @ model_design[rows]
                score += model_design[rows].T @ weight @ local_values
                quadratic += local_values @ weight @ local_values
                log_determinant += np.linalg.slogdet(variance)[1]
            coefficient = np.linalg.solve(precision, score)
            return log_determinant + quadratic - score @ coefficient

        self.assertAlmostEqual(
            result["null_objective"],
            dense_profile_objective(
                np.ones((len(values), 1)),
                result["null_cluster_variance"],
                result["null_residual_variance"],
            ),
            places=8,
        )
        self.assertAlmostEqual(
            result["alternative_objective"],
            dense_profile_objective(
                np.column_stack((np.ones(len(values)), treatment)),
                result["alternative_cluster_variance"],
                result["alternative_residual_variance"],
            ),
            places=8,
        )
        dense = differential.gaussian_mixed_model_lrt(
            values,
            covariances,
            np.column_stack((np.ones(len(values)), treatment)),
            tested_columns=[1],
            clusters=clusters,
            effect_prior_scales=(0.25, 0.5, 1.0),
            random_intercept_solver="dense",
        )
        for key in (
            "null_objective",
            "alternative_objective",
            "statistic",
            "mixture_log_bayes_factor",
        ):
            self.assertAlmostEqual(result[key], dense[key], places=8)

    def test_gaussian_score_bf_penalizes_zero_effect(self):
        subjects = 12
        clusters = np.repeat(np.arange(subjects), 2)
        treatment = np.tile([0.0, 1.0], subjects)
        values = np.zeros((len(treatment), 2))
        covariances = np.repeat(
            np.eye(2)[None, :, :] * 0.1, len(values), axis=0
        )
        result = differential.gaussian_mixed_model_lrt(
            values,
            covariances,
            np.column_stack((np.ones(len(values)), treatment)),
            tested_columns=[1],
            clusters=clusters,
            effect_prior_scales=(0.25,),
        )
        self.assertAlmostEqual(result["statistic"], 0.0)
        self.assertLess(result["mixture_log_bayes_factor"], 0.0)
        self.assertEqual(result["log_bayes_factors"].shape, (1,))

    def test_gaussian_score_bf_is_invariant_to_reference_level(self):
        clusters = np.repeat(np.arange(12), 3)
        levels = np.tile(np.arange(3), 12)
        values = np.tile(np.array([0.0, 0.3, -0.2]), 12)[:, None]
        covariances = np.repeat(
            np.array([[[0.1]]]), len(values), axis=0
        )
        prior_covariance = np.eye(2) + np.ones((2, 2))
        first = differential.gaussian_mixed_model_lrt(
            values,
            covariances,
            np.column_stack((
                np.ones(len(values)), levels == 1, levels == 2
            )),
            tested_columns=[1, 2],
            clusters=clusters,
            effect_prior_scales=(0.25,),
            effect_prior_covariance=prior_covariance,
        )
        second = differential.gaussian_mixed_model_lrt(
            values,
            covariances,
            np.column_stack((
                np.ones(len(values)), levels == 0, levels == 2
            )),
            tested_columns=[1, 2],
            clusters=clusters,
            effect_prior_scales=(0.25,),
            effect_prior_covariance=prior_covariance,
        )
        self.assertAlmostEqual(
            first["mixture_log_bayes_factor"],
            second["mixture_log_bayes_factor"],
            places=8,
        )

    def test_paired_cell_type_test_detects_effect(self):
        rng = np.random.default_rng(9)
        records = []
        for mouse in range(10):
            condition = "case" if mouse >= 5 else "control"
            for cell_type, mean in (("A", 0.0), ("B", 1.5)):
                records.append(
                    (
                        np.array([mean + rng.normal(0, 0.05)]),
                        np.array([[0.01]]),
                        cell_type,
                        condition,
                        f"mouse_{mouse}",
                    )
                )
        inputs = {("block", "conditional"): records}
        with tempfile.TemporaryDirectory() as directory:
            args = type(
                "Args",
                (),
                {
                    "min_celltype_mice": 3,
                    "permutations": 5,
                    "output_dir": Path(directory),
                },
            )()
            summary = (
                run_differential_splicing.cell_type_differential_tests(
                    args,
                    inputs,
                    np.random.default_rng(10),
                )
            )
            table = np.genfromtxt(
                Path(directory) / "differential_cell_type.tsv",
                delimiter="\t",
                names=True,
                dtype=None,
                encoding="utf-8",
            )
        self.assertEqual(summary["cell_type_tests"], 1)
        self.assertEqual(summary["cell_type_fdr_blocks"], 1)
        self.assertLess(float(table["p_value"]), 0.01)


if __name__ == "__main__":
    unittest.main()
