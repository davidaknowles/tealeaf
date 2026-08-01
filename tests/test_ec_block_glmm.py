"""Tests for block-specific fixed effects in EC-count GLMMs."""

import numpy as np

from tealeaf.sc import ec_block_glmm, ec_glmm, ec_glmm_full


def test_path_contrasts_tie_isoforms_on_the_same_path():
    contrasts = ec_block_glmm.path_contrast_matrix([0, 0, 1, -1])
    assert contrasts.shape == (3, 1)
    np.testing.assert_allclose(contrasts[0], contrasts[1])
    assert not np.allclose(contrasts[0], contrasts[2])


def test_block_tensors_are_nested_with_expected_df():
    nuisance = np.array([[1, 0], [0, 1], [1, 0]], dtype=float)
    tested = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    null, alternative, degrees = ec_block_glmm.block_fixed_effect_tensors(
        nuisance, tested, [0, 0, 1, 2]
    )
    assert null.shape == (3, 3, 6)
    assert alternative.shape == (3, 3, 10)
    assert degrees == 4
    np.testing.assert_array_equal(alternative[:, :, :6], null)


def test_subject_permutation_preserves_each_subject_labels():
    values = np.array([0, 1, 2, 0, 1, 2])
    subjects = np.array(["a", "a", "a", "b", "b", "b"])
    result = ec_block_glmm.permute_within_subject(
        values, subjects, np.random.default_rng(4)
    )
    for subject in np.unique(subjects):
        positions = subjects == subject
        np.testing.assert_array_equal(
            np.sort(result[positions]), np.sort(values[positions])
        )


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
