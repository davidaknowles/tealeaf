import numpy as np
import pytest
from scipy.special import logit

from tealeaf.sc.local_read_reporting import pooled_local_marker_effect


def test_expression_weights_can_reverse_pooled_usage_despite_identical_positive_subject_log_odds():
    counts = np.array([[[[1, 9], [200, 800]]], [[[800, 200], [9, 1]]]], dtype=np.int64)
    proportions = counts[..., 0] / counts.sum(axis=-1)
    effects = logit(proportions[:, 0, 1]) - logit(proportions[:, 0, 0])
    assert np.allclose(effects, effects[0]) and (effects > 0).all()
    result = pooled_local_marker_effect(counts)
    assert np.isclose(result['pooled_marker_effect'], 209 / 1010 - 801 / 1010)
    assert result['pooled_marker_effect'] < 0 and result['pooled_marker_n_primers'] == 1


def test_primers_are_standardized_not_confounded_by_between_primer_depth():
    # Within each primer the effect is positive; naive pooling reverses it.
    counts = np.array([[[[800, 200], [9, 1]], [[1, 9], [200, 800]]]], dtype=np.int64)
    result = pooled_local_marker_effect(counts)
    assert np.isclose(result['pooled_marker_effect'], .1)
    assert result['pooled_marker_n_primers'] == 2


def test_absent_strata_and_subject_reordering_do_not_invent_information():
    counts = np.array([[[[1, 9], [2, 8]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[0, 0], [3, 7]]]], dtype=np.int64)
    result = pooled_local_marker_effect(counts)
    assert np.isclose(result['pooled_marker_effect'], .1)
    assert result['pooled_marker_n_primers'] == 1 and result['pooled_marker_primer_effects'][1] is None
    assert pooled_local_marker_effect(counts[::-1]) == result
    missing = pooled_local_marker_effect(np.zeros_like(counts))
    assert np.isnan(missing['pooled_marker_effect']) and missing['pooled_marker_n_primers'] == 0


@pytest.mark.parametrize('bad', [np.nan, -1, .25])
def test_reporting_rejects_invalid_original_counts(bad):
    counts = np.ones((4, 2, 2, 2))
    counts[0, 0, 0, 0] = bad
    with pytest.raises(ValueError):
        pooled_local_marker_effect(counts)
