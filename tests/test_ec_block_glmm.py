"""Tests for block-specific fixed effects in EC-count GLMMs."""

import numpy as np
import pandas as pd

from extra_scripts.merge_ec_block_glmm import calibrate
from extra_scripts.run_ec_block_glmm import (
    covered_celltype_design,
    covered_condition_designs,
    modeled_gene_umis,
)
from tealeaf.sc import ec_block_glmm, ec_glmm, ec_glmm_full


def test_path_contrasts_tie_isoforms_on_the_same_path():
    contrasts = ec_block_glmm.path_contrast_matrix([0, 0, 1, -1])
    assert contrasts.shape == (3, 1)
    np.testing.assert_allclose(contrasts[0], contrasts[1])
    assert not np.allclose(contrasts[0], contrasts[2])


def test_path_and_nuisance_bases_span_centered_isoform_logits():
    path, nuisance = ec_block_glmm.block_effect_bases([0, 0, 1, -1])
    assert path.shape == (3, 1)
    assert nuisance.shape == (3, 2)
    assert np.linalg.matrix_rank(np.column_stack((nuisance, path))) == 3

    def full_centered(reference_logits):
        full = np.vstack((reference_logits, np.zeros(reference_logits.shape[1])))
        return full - full.mean(axis=0, keepdims=True)

    np.testing.assert_allclose(
        full_centered(path).T @ full_centered(nuisance), 0.0, atol=1e-12
    )


def test_block_tensors_are_nested_with_expected_df():
    nuisance = np.array([[1, 0], [0, 1], [1, 0]], dtype=float)
    tested = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    null, alternative, degrees = ec_block_glmm.block_fixed_effect_tensors(
        nuisance, tested, [0, 0, 1, 2]
    )
    assert null.shape == (3, 3, 8)
    assert alternative.shape == (3, 3, 12)
    assert degrees == 4
    np.testing.assert_array_equal(alternative[:, :, :8], null)


def test_full_vi_accepts_tensor_fixed_effects_and_warm_start():
    counts = (np.array([[8, 2], [3, 7], [7, 3], [2, 8]], dtype=float),)
    compatibility = (np.eye(2),)
    design = np.column_stack((np.ones(4), [0, 1, 0, 1]))
    clusters = np.array([0, 0, 1, 1])
    tensor = ec_block_glmm.unrestricted_tensor(design[:, :1], 1)
    null_data = ec_glmm.ECGLMMData(
        counts, compatibility, np.ones((4, 1)), clusters,
        fixed_effect_tensor=tensor,
    )
    null = ec_glmm_full.fit_variational(null_data, max_iter=3)
    expanded = np.concatenate(
        (tensor, ec_block_glmm.unrestricted_tensor(design[:, 1:], 1)), axis=2
    )
    alternative_data = ec_glmm.ECGLMMData(
        counts, compatibility, np.ones((4, 1)), clusters,
        fixed_effect_tensor=expanded,
    )
    initial = ec_glmm_full.fixed_effect_warm_start(null, 2)
    alternative = ec_glmm_full.fit_variational(
        alternative_data, initial=initial, max_iter=3
    )
    assert np.isfinite(alternative["objective"])
    assert alternative["fixed_effect_count"] == 2


def test_null_simulation_preserves_primer_totals():
    counts = (
        np.array([[8, 2], [3, 7]], dtype=float),
        np.array([[4, 6], [9, 1]], dtype=float),
    )
    tensor = np.ones((2, 1, 1), dtype=float)
    data = ec_glmm.ECGLMMData(
        counts,
        (np.eye(2), np.eye(2)),
        np.ones((2, 1)),
        np.array(["a", "b"]),
        fixed_effect_tensor=tensor,
    )
    fit = {
        "coefficients": np.array([0.3]),
        "random_effect_sd": 0.2,
        "observation_noise_sd": 0.1,
        "concentration": 15.0,
    }
    simulated = ec_block_glmm.simulate_null_counts(
        data,
        fit,
        np.random.default_rng(3),
        family="dirichlet_multinomial",
        observation_noise=True,
    )
    for original, generated in zip(counts, simulated):
        np.testing.assert_array_equal(
            original.sum(axis=1), generated.sum(axis=1)
        )


def test_null_simulation_accepts_primer_without_modeled_ecs():
    counts = (np.array([[8, 2], [3, 7]], dtype=float), np.empty((2, 0)))
    data = ec_glmm.ECGLMMData(
        counts,
        (np.eye(2), np.empty((0, 2))),
        np.ones((2, 1)),
        np.array(["a", "b"]),
        fixed_effect_tensor=np.ones((2, 1, 1)),
    )
    fit = {
        "coefficients": np.array([0.3]),
        "random_effect_sd": 0.2,
        "observation_noise_sd": 0.0,
        "concentration": np.inf,
    }
    simulated = ec_block_glmm.simulate_null_counts(
        data,
        fit,
        np.random.default_rng(3),
        family="multinomial",
    )
    assert simulated[1].shape == (2, 0)


def test_bootstrap_calibration_leaves_out_tested_block():
    table = pd.DataFrame({
        "test_id": ["b1", "b2"],
        "block_id": ["b1", "b2"],
        "method": ["multinomial_full"] * 2,
        "df_stratum": ["1"] * 2,
        "depth_bin": [0] * 2,
        "statistic": [2.5, 0.5],
    })
    null = pd.DataFrame({
        "test_id": ["b1", "b1", "b2", "b2"],
        "block_id": ["b1", "b1", "b2", "b2"],
        "method": ["multinomial_full"] * 4,
        "replicate": [0, 1, 0, 1],
        "depth_bin": [0] * 4,
        "statistic": [0.0, 1.0, 2.0, 3.0],
        "converged": [True] * 4,
    })
    calibrated, calibrated_null = calibrate(table, null)
    np.testing.assert_allclose(calibrated["p_value"], [2 / 3, 2 / 3])
    assert len(calibrated_null) == 4


def test_covered_celltype_design_is_gene_specific():
    metadata = pd.DataFrame({
        "cell_type": ["a"] * 3 + ["b"] * 3 + ["c"] * 3,
        "condition": ["x", "x", "y"] * 3,
        "mouse": ["m1", "m2", "m3"] * 3,
    })
    result = covered_celltype_design(
        metadata,
        [25, 30, 40, 50, 30, 35, 100, 100, 0],
        min_gene_umis=25,
        min_samples=6,
        min_cell_types=2,
        min_celltype_mice=3,
    )
    rows, local, nuisance, labels, _, cell_types = result
    np.testing.assert_array_equal(rows, np.arange(6))
    assert cell_types == ["a", "b"]
    assert nuisance.shape == (6, 2)
    np.testing.assert_array_equal(labels, [0, 0, 0, 1, 1, 1])
    assert set(local["cell_type"]) == {"a", "b"}


def test_condition_designs_are_built_within_cell_type():
    metadata = pd.DataFrame({
        "cell_type": ["a"] * 6 + ["b"] * 5,
        "condition": ["x"] * 3 + ["y"] * 3 + ["x"] * 3 + ["y"] * 2,
        "mouse": ["m1", "m2", "m3"] * 3 + ["m4", "m5"],
    })
    results = covered_condition_designs(
        metadata,
        np.full(len(metadata), 25),
        min_gene_umis=10,
        min_samples=0,
        min_conditions=2,
        min_condition_mice=3,
    )
    assert len(results) == 1
    coverage, cell_type = results[0]
    rows, local, nuisance, labels, _, conditions = coverage
    assert cell_type == "a"
    assert conditions == ["x", "y"]
    np.testing.assert_array_equal(rows, np.arange(6))
    np.testing.assert_array_equal(labels, [0, 0, 0, 1, 1, 1])
    np.testing.assert_array_equal(nuisance, np.ones((6, 1)))
    assert set(local["cell_type"]) == {"a"}


def test_modeled_gene_umis_exclude_unsupported_primer_ecs():
    counts = (
        np.array([[5, 7], [11, 13]]),
        np.array([[17, 19], [23, 29]]),
    )
    designs = (
        np.array([[1, 0], [0, 0]]),
        np.array([[0, 0], [0, 1]]),
    )
    np.testing.assert_array_equal(
        modeled_gene_umis(counts, designs, [0, 1], [0, 1]),
        [24, 40],
    )
