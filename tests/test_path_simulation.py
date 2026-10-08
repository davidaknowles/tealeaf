import numpy as np
import pytest

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_simulation import simulate_counts, resample_counts


def test_latent_path_draw_is_shared_across_technical_observations_and_primers():
    class DeterministicRng:
        calls = 0

        def normal(self, scale, size):
            return np.zeros(size)

        def dirichlet(self, parameters):
            self.calls += 1
            return np.asarray([.2, .8] if self.calls % 2 else [.8, .2])

        def multinomial(self, total, probabilities):
            return total * probabilities

    subjects = np.repeat([0, 1], 4)
    labels = np.tile([0, 0, 1, 1], 2)
    counts = np.tile([20., 30., 50.], (8, 1))
    data = ECGLMMData((counts, counts), (np.eye(3), np.diag([2., 1., 1.])), np.ones((8, 1)), subjects)
    rng = DeterministicRng()
    result = simulate_counts(data, [.2, .3, .5], subjects, rng, subject_scale=0, labels=labels, path_index=[0, 1, -1], residual_concentration=20)
    assert rng.calls == 4
    assert np.allclose(result.counts[0][labels == 0], [10., 40., 50.])
    assert np.allclose(result.counts[0][labels == 1], [40., 10., 50.])
    assert np.allclose(result.counts[1][labels == 0], np.asarray([20., 40., 50.]) * (100 / 110))
    assert np.allclose(result.counts[1][labels == 1], np.asarray([80., 10., 50.]) * (100 / 140))
    assert all(np.allclose(matrix.sum(axis=1), 100) for matrix in result.counts)


def test_finite_dispersion_has_equal_expected_path_usage_across_labels():
    subjects = np.zeros(2, int)
    labels = np.arange(2)
    data = ECGLMMData((np.tile([40., 60.], (2, 1)),), (np.eye(2),), np.ones((2, 1)), subjects)
    rng = np.random.default_rng(7)
    draws = np.asarray([simulate_counts(data, [.4, .6], subjects, rng, subject_scale=0, labels=labels, path_index=[0, 1], residual_concentration=20).counts[0][:, 0] / 100 for _ in range(1000)])
    assert np.allclose(draws.mean(axis=0), .4, atol=.025)
    # Includes biological Dirichlet variance, not just read-count variance.
    assert np.all(draws.var(axis=0) > .007)


def test_residual_variation_requires_valid_paths_and_labels():
    data = ECGLMMData((np.ones((4, 2)),), (np.eye(2),), np.ones((4, 1)), np.arange(4))
    with pytest.raises(ValueError, match="aligned labels"):
        simulate_counts(data, [.4, .6], np.arange(4), np.random.default_rng(0), residual_concentration=20)
    with pytest.raises(ValueError, match="positive and finite"):
        simulate_counts(data, [.4, .6], np.arange(4), np.random.default_rng(0), labels=np.arange(4), path_index=[0, 1], residual_concentration=0)


def test_truth_export_preserves_default_random_draws_and_fixed_composition():
    subjects, labels = np.repeat(np.arange(4), 2), np.tile([0, 1], 4)
    data = ECGLMMData((np.tile([40., 60.], (8, 1)), np.tile([50., 70.], (8, 1))), (np.eye(2), np.diag([2., 1.])), np.ones((8, 1)), subjects)
    options = {"labels": labels, "path_index": [0, 1], "residual_concentration": 20}
    first = simulate_counts(data, [.4, .6], subjects, np.random.default_rng(87), **options)
    second, details = simulate_counts(data, [.4, .6], subjects, np.random.default_rng(87), return_details=True, **options)
    assert all(np.array_equal(a, b) for a, b in zip(first.counts, second.counts))
    assert details["observation_weights"].shape == (8, 2)
    assert np.allclose(details["subject_weights"].sum(axis=1), 1)
    repeated = resample_counts(data, details["observation_weights"], np.random.default_rng(5))
    assert all(np.array_equal(a.sum(axis=1), b.sum(axis=1)) for a, b in zip(repeated.counts, data.counts))
    assert not np.array_equal(repeated.counts[0], first.counts[0])


def test_fixed_composition_resampling_rejects_bad_weights():
    data = ECGLMMData((np.ones((4, 2)),), (np.eye(2),), np.ones((4, 1)), np.arange(4))
    for weights in (np.ones((3, 2)), np.zeros((4, 2)), -np.ones((4, 2))):
        with pytest.raises(ValueError, match="observation-by-transcript"):
            resample_counts(data, weights, np.random.default_rng(0))
def test_independent_binary_blocks_preserve_null_A_and_primer_totals():
    from tealeaf.sc.path_simulation import simulate_independent_binary_blocks

    data, truth = simulate_independent_binary_blocks(np.random.default_rng(5), n_subjects=6, b_effect=2.2)
    assert data.n_isoforms == 4
    assert truth["true_delta"] == 0.
    np.testing.assert_array_equal(truth["a_usage"][::2], truth["a_usage"][1::2])
    assert np.all(truth["b_usage"][::2] < .11)
    assert np.all(truth["b_usage"][1::2] > .89)
    np.testing.assert_allclose(truth["observation_weights"][:, :2].sum(axis=1), truth["a_usage"])
    for primer, counts in enumerate(data.counts):
        assert counts.shape == (12, 7)
        np.testing.assert_array_equal(counts.sum(axis=1), 100 * (primer + 1))
    again, _ = simulate_independent_binary_blocks(np.random.default_rng(5), n_subjects=6, b_effect=2.2)
    for first, second in zip(data.counts, again.counts):
        np.testing.assert_array_equal(first, second)
def test_within_path_type_tilts_preserve_local_null_and_outside_mass():
    baseline = np.array([.15, .25, .4, .2])
    subjects, labels = np.repeat(np.arange(6), 2), np.tile([0, 1], 6)
    counts = np.tile([3, 5, 8, 4], (12, 1))
    base = ECGLMMData((counts,), (np.eye(4),), np.ones((12, 1)), subjects)
    _, truth = simulate_counts(base, baseline, subjects, np.random.default_rng(31), labels=labels, path_index=np.array([0, 0, 1, -1]), within_path_type_scale=2., return_details=True)
    weights = truth["observation_weights"]
    np.testing.assert_allclose(weights[::2, :2].sum(axis=1), weights[1::2, :2].sum(axis=1))
    np.testing.assert_allclose(weights[::2, 2:], weights[1::2, 2:])
    assert not np.allclose(weights[::2, 0] / weights[::2, :2].sum(axis=1), weights[1::2, 0] / weights[1::2, :2].sum(axis=1))


@pytest.mark.parametrize("residual", [None, 20.])
def test_event_mass_type_changes_preserve_conditional_path_truth(residual):
    subjects = np.repeat(np.arange(6), 4)
    labels = np.tile([0, 0, 1, 1], 6)
    counts = np.tile([3., 5., 8., 4.], (24, 1))
    base = ECGLMMData((counts, counts * 2), (np.eye(4), np.diag([2., 1., 3., 1.])), np.ones((24, 1)), subjects)
    options = dict(labels=labels, path_index=[0, 0, 1, -1], residual_concentration=residual, within_path_type_scale=1., return_details=True)
    _, original = simulate_counts(base, [.15, .25, .4, .2], subjects, np.random.default_rng(31), **options)
    generated, shifted = simulate_counts(base, [.15, .25, .4, .2], subjects, np.random.default_rng(31), event_mass_type_scale=2., **options)
    before, after = original["observation_weights"], shifted["observation_weights"]
    np.testing.assert_allclose(after[:, :3] / after[:, :3].sum(axis=1, keepdims=True), before[:, :3] / before[:, :3].sum(axis=1, keepdims=True))
    np.testing.assert_allclose(after.sum(axis=1), 1.)
    assert not np.allclose(after[:, :3].sum(axis=1), before[:, :3].sum(axis=1))
    np.testing.assert_array_equal(shifted["event_mass_type_levels"], [0, 1])
    np.testing.assert_allclose(shifted["event_mass_type_tilts"].sum(), 0., atol=1e-15)
    np.testing.assert_array_equal(after[::4], after[1::4])
    np.testing.assert_array_equal(after[2::4], after[3::4])
    for primer, matrix in enumerate(generated.counts):
        np.testing.assert_array_equal(matrix.sum(axis=1), 20. * (primer + 1))


@pytest.mark.parametrize("scale,paths", [(None, [0, 1, -1]), (0., [0, 1, -1]), (2., [0, 0, 1])])
def test_disabled_or_inapplicable_event_mass_changes_preserve_random_stream(scale, paths):
    subjects = np.repeat(np.arange(4), 2)
    counts = np.tile([2., 3., 5.], (8, 1))
    base = ECGLMMData((counts,), (np.eye(3),), np.ones((8, 1)), subjects)
    options = dict(labels=np.tile([0, 1], 4), path_index=paths, return_details=True)
    first_rng, second_rng = np.random.default_rng(813), np.random.default_rng(813)
    first, original = simulate_counts(base, [.2, .3, .5], subjects, first_rng, **options)
    second, shifted = simulate_counts(base, [.2, .3, .5], subjects, second_rng, event_mass_type_scale=scale, **options)
    np.testing.assert_array_equal(first.counts[0], second.counts[0])
    np.testing.assert_array_equal(original["observation_weights"], shifted["observation_weights"])
    assert first_rng.bit_generator.state == second_rng.bit_generator.state


@pytest.mark.parametrize("scale", [-1., np.inf, np.nan])
def test_event_mass_scale_must_be_finite_and_nonnegative(scale):
    base = ECGLMMData((np.ones((4, 3)),), (np.eye(3),), np.ones((4, 1)), np.arange(4))
    with pytest.raises(ValueError, match="finite nonnegative"):
        simulate_counts(base, [.2, .3, .5], np.arange(4), np.random.default_rng(3), event_mass_type_scale=scale)


def test_event_mass_changes_require_aligned_labels_and_paths():
    base = ECGLMMData((np.ones((4, 3)),), (np.eye(3),), np.ones((4, 1)), np.arange(4))
    for options in ({}, dict(labels=[0, 1], path_index=[0, 1, -1]), dict(labels=[0, 1, 0, 1], path_index=[0, 1])):
        with pytest.raises(ValueError, match="aligned labels and paths"):
            simulate_counts(base, [.2, .3, .5], np.arange(4), np.random.default_rng(3), event_mass_type_scale=1., **options)
