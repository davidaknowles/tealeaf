import numpy as np
import pytest

from tealeaf.sc.event_paths import binary_event_information, mixture_event_information


def test_explicit_mixture_information_handles_unequal_depth_and_empty_primer():
    result = mixture_event_information((np.eye(2), np.zeros((3, 2))), [0, 1], [.5, .5], [[100, 200], [0, 0]])
    assert result["free_transcript_information"] == pytest.approx(2 * .5 * .5 * 100 * 200 / 300)
    assert result["fixed_mixture_information"] == pytest.approx(result["free_transcript_information"])


def test_identifiable_categorical_event_keeps_information_without_fixed_shares():
    result = binary_event_information((np.eye(3),), [0, 0, 1])
    assert result["free_transcript_information"] > 0
    assert result["free_transcript_information"] == pytest.approx(result["fixed_mixture_information"])


def test_fixed_shares_can_manufacture_identifiability_of_aliased_event():
    mapping = np.array([[1., 0, 0], [0, 1, 1]])
    result = binary_event_information((mapping,), [0, 0, 1])
    assert result["free_transcript_information"] < 1e-10
    assert result["fixed_mixture_information"] > 100


def test_both_assumptions_reject_information_in_identical_event_columns():
    result = binary_event_information((np.ones((2, 2)),), [0, 1])
    assert result["free_transcript_information"] < 1e-10
    assert result["fixed_mixture_information"] < 1e-10


def test_information_audit_requires_positive_primer_support_and_interior_paths():
    with pytest.raises(ValueError, match="mass"):
        binary_event_information((np.zeros((2, 2)),), [0, 1])
    with pytest.raises(ValueError, match="interior"):
        binary_event_information((np.eye(2),), [0, 1], inclusion=0)
