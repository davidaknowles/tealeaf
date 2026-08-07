"""Tests for block-specific fixed effects in EC-count GLMMs."""

from types import SimpleNamespace

import numpy as np
import pandas as pd

from extra_scripts.merge_ec_block_glmm import (
    calibrate,
    calibrate_gpd_tail,
    fit_gpd_tail,
    gpd_tail_p_values,
)
from extra_scripts.run_ec_block_glmm import (
    covered_celltype_design,
    covered_celltype_pairwise_designs,
    covered_condition_designs,
    deduplicate_supported_partitions,
    fit_method,
    modeled_gene_umis,
    joint_gene_candidates,
    partition_candidates,
    permute_paired_labels,
    reparameterize_fixed_effects,
)
from extra_scripts.run_celltype_compositional_splicing import (
    paired_celltype_design,
)
from extra_scripts.run_differential_splicing import block_mapping
from tealeaf.sc import ec_block_glmm, ec_glmm, ec_glmm_full


def test_laplace_dirichlet_multinomial_dispatches_family(monkeypatch):
    captured = {}

    def fake_fit(data, **kwargs):
        captured.update(kwargs)
        return {"data": data}

    monkeypatch.setattr(ec_glmm, "fit_laplace", fake_fit)
    result = fit_method(
        "laplace_dirichlet_multinomial",
        "sentinel",
        SimpleNamespace(max_iter=17, vi_samples=3, seed=1),
    )
    assert result["data"] == "sentinel"
    assert captured["family"] == "dirichlet_multinomial"
    assert captured["observation_noise"] is False
    assert captured["random_slopes"] is False
    assert captured["max_iter"] == 17


def test_celltype_design_can_omit_subject_fixed_effects():
    records = [
        {"cluster": cell_type, "condition": "A", "mouse": mouse}
        for mouse in ("m1", "m2", "m3")
        for cell_type in ("x", "y")
    ]
    fixed = paired_celltype_design(records, 3, "fixed")
    marginal = paired_celltype_design(records, 3, "none")
    assert fixed[-2].shape == (6, 3)
    assert fixed[-1].shape == (6, 4)
    assert marginal[-2].shape == (6, 1)
    assert marginal[-1].shape == (6, 2)
    np.testing.assert_array_equal(marginal[-2], 1.0)


def test_block_mapping_ignores_ensembl_version_suffixes():
    block = SimpleNamespace(
        transcripts=("ENST1.4", "ENST2.1"),
        path_index=(0, 1),
        path_signatures=(((),), ((),)),
    )
    mapping, _ = block_mapping(block, np.array(["ENST1", "ENST2"]))
    np.testing.assert_array_equal(mapping, [0, 1])


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


def test_canonical_full_alternative_has_same_block_model_span():
    nuisance = np.array([[1, 0], [0, 1], [1, 0]], dtype=float)
    tested = np.array([[0], [1], [0]], dtype=float)
    null, block_alternative, _ = ec_block_glmm.block_fixed_effect_tensors(
        nuisance, tested, [0, 0, 1, 2]
    )
    canonical = ec_block_glmm.full_fixed_effect_tensor(nuisance, tested, 3)
    block_flat = block_alternative.reshape(-1, block_alternative.shape[2])
    canonical_flat = canonical.reshape(-1, canonical.shape[2])
    assert np.linalg.matrix_rank(block_flat) == np.linalg.matrix_rank(canonical_flat)
    assert np.linalg.matrix_rank(np.column_stack((block_flat, canonical_flat))) == (
        np.linalg.matrix_rank(canonical_flat)
    )

    parameters = np.r_[np.arange(null.shape[2], dtype=float), -0.7]
    initial = reparameterize_fixed_effects(
        {"parameters": parameters}, null, canonical
    )
    np.testing.assert_allclose(
        canonical_flat @ initial[: canonical.shape[2]],
        null.reshape(-1, null.shape[2]) @ parameters[:-1],
    )
    assert initial[-1] == parameters[-1]


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


def test_pooled_weights_and_path_collapse_reduce_nuisance_dimension():
    mappings = (
        np.array([[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]),
    )
    counts = (np.array([[80, 20, 30, 10], [80, 20, 30, 10]]),)
    weights = ec_block_glmm.pooled_isoform_weights(counts, mappings)
    assert weights[0] > weights[1]
    collapsed, paths, projection = ec_block_glmm.collapse_within_paths(
        mappings, weights, np.array([0, 0, 1, -1])
    )
    assert collapsed[0].shape == (4, 3)
    np.testing.assert_array_equal(paths, [0, 1, -1])
    assert projection[0, 0] > projection[1, 0]
    np.testing.assert_allclose(projection.sum(axis=0), 1.0)


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


def test_pairwise_celltype_designs_require_gene_covered_shared_subjects():
    metadata = pd.DataFrame({
        "cell_type": ["a", "a", "a", "b", "b", "b", "c", "c"],
        "condition": ["x", "x", "y", "x", "x", "y", "x", "y"],
        "mouse": ["m1", "m2", "m3", "m1", "m2", "m4", "m1", "m3"],
    })
    results = covered_celltype_pairwise_designs(
        metadata,
        np.full(len(metadata), 25),
        min_gene_umis=10,
        min_samples=0,
        min_celltype_mice=2,
    )
    assert len(results) == 2
    by_levels = {coverage[-1]: coverage for coverage, _ in results}
    paired = by_levels[("a", "b")]
    rows, local, nuisance, labels, subjects, levels = paired
    np.testing.assert_array_equal(rows, [0, 1, 3, 4])
    np.testing.assert_array_equal(labels, [0, 0, 1, 1])
    assert levels == ("a", "b")
    assert set(subjects) == {"m1", "m2"}
    assert nuisance.shape == (4, 1)
    assert set(local.cell_type) == {"a", "b"}


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


def test_supported_partition_deduplication_ignores_path_labels():
    common = ("gene", 0, np.array([1, 2, 3]), None)
    suffix = ((), np.array([0, 1]), "ct", ("a", "b"))
    first = (
        "b1", "b1", *common[:2], common[2], np.array([0, 0, 1]), *suffix
    )
    duplicate = (
        "b2", "b2", *common[:2], common[2], np.array([2, 2, 5]), *suffix
    )
    distinct = (
        "b3", "b3", *common[:2], common[2], np.array([0, 1, 1]), *suffix
    )
    assert [row[0] for row in deduplicate_supported_partitions(
        [first, duplicate, distinct]
    )] == ["b1", "b3"]


def test_candidate_partition_keeps_shared_alternatives_together():
    def candidate(test, gene, rows):
        return (
            test, test, gene, gene, np.array([0, 1]), np.array([0, 1]),
            (), np.array(rows), None, ("a", "b"),
        )

    candidates = [
        candidate("a1", 1, [0, 1]),
        candidate("a2", 1, [0, 1]),
        candidate("a3", 1, [0, 1]),
        candidate("b1", 2, [2, 3]),
        candidate("b2", 2, [2, 3]),
        candidate("c1", 3, [4, 5]),
    ]
    shards = partition_candidates(candidates, 2)
    locations = {
        row[0]: index for index, shard in enumerate(shards) for row in shard
    }
    assert locations["a1"] == locations["a2"] == locations["a3"]
    assert locations["b1"] == locations["b2"]
    assert sorted(map(len, shards)) == [3, 3]


def test_joint_gene_candidates_keep_one_unrestricted_test_per_design():
    common = (
        "g1.1", 1, np.array([4, 5, 6]), None, np.array([0, 1]), None,
        ("a", "b"),
    )
    candidates = [
        ("b1", "b1", common[0], common[1], common[2], np.array([0, 0, 1]),
         common[3], common[4], common[5], common[6]),
        ("b2", "b2", common[0], common[1], common[2], np.array([0, 1, 1]),
         common[3], common[4], common[5], common[6]),
    ]
    result = joint_gene_candidates(candidates)
    assert len(result) == 1
    assert result[0][0] == "g1.1|joint|a|b"
    assert result[0][1] == "g1.1:JOINT"
    np.testing.assert_array_equal(result[0][5], [0, 1, 2])


def test_pairwise_null_permutation_swaps_only_within_subject():
    metadata = pd.DataFrame({
        "mouse": np.repeat([f"m{i}" for i in range(12)], 2),
        "cell_type": ["a", "b"] * 12,
    })
    labels = np.tile([0, 1], 12)
    permuted = permute_paired_labels(metadata, labels, ("a", "b"), 7)
    np.testing.assert_array_equal(
        permuted.reshape(-1, 2).sum(axis=1), np.ones(12, dtype=int)
    )
    assert np.any(permuted != labels)
    assert np.any(permuted == labels)
    np.testing.assert_array_equal(
        permuted,
        permute_paired_labels(metadata, labels, ("a", "b"), 7),
    )


def test_gpd_tail_resolves_beyond_empirical_bootstrap_floor():
    values = np.random.default_rng(4).exponential(size=1000)
    model = fit_gpd_tail(values, 0.9)
    p_value = gpd_tail_p_values([values.max() + 5], model)[0]
    assert 0 < p_value < 1 / (1 + len(values))


def test_gpd_tail_calibration_cross_fits_whole_tests():
    rng = np.random.default_rng(8)
    test_ids = [f"b{index}" for index in range(60)]
    table = pd.DataFrame({
        "test_id": test_ids,
        "block_id": test_ids,
        "method": "multinomial_full",
        "df_stratum": "1",
        "depth_bin": 0,
        "statistic": rng.exponential(size=len(test_ids)),
    })
    null = pd.DataFrame({
        "test_id": np.repeat(test_ids, 5),
        "block_id": np.repeat(test_ids, 5),
        "method": "multinomial_full",
        "replicate": np.tile(np.arange(5), len(test_ids)),
        "depth_bin": 0,
        "statistic": rng.exponential(size=5 * len(test_ids)),
        "converged": True,
    })
    calibrated, calibrated_null = calibrate_gpd_tail(
        table, null, quantile=0.9, folds=5
    )
    assert calibrated["p_value"].between(0, 1).all()
    assert calibrated_null["p_value"].between(0, 1).all()
    assert len(calibrated_null) == len(null)
