import numpy as np
import pytest

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_simulation import simulate_counts


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
