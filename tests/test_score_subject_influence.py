import numpy as np
import pytest

from tealeaf.sc.path_score_mixed import PathScoreComponents, binary_score_components_to_proportions, binary_score_subject_influence, efficient_shared_path_score


def example_components(information):
    information = np.asarray(information, float)
    size = len(information)
    return PathScoreComponents((information * .5)[:, None], information[:, None, None], np.ones((size, 1, 1)), np.array([f"s{index}" for index in range(size)]), (0, 1), [], [[] for _ in range(size)], "ilr", np.maximum(information, 1.)[:, None, None])


def test_equal_subject_information_gives_equal_precision_without_mutation():
    components = example_components([4.] * 6)
    before = [value.copy() for value in (components.scores, components.information, components.biological_shapes, components.reference_information)]
    records, summary = binary_score_subject_influence(components)
    np.testing.assert_allclose([row["precision_share"] for row in records], 1 / 6)
    assert summary["effective_weighted_subjects"] == pytest.approx(6.)
    assert summary["fitted_mean"] == pytest.approx(.5)
    assert summary["stable_fixed_variance_p_value"] == pytest.approx(summary["fitted_p_value"], rel=1e-10)
    for original, value in zip(before, (components.scores, components.information, components.biological_shapes, components.reference_information)):
        np.testing.assert_array_equal(original, value)


def test_one_high_information_subject_is_not_counted_as_five_effective_subjects():
    records, summary = binary_score_subject_influence(example_components([1e6, 1., 1., 1., 1.]))
    assert summary["maximum_precision_share"] > .9999
    assert summary["effective_weighted_subjects"] < 1.001
    assert summary["dominant_subject"] == "s0"
    assert len(records) == 5 and all(row["retained"] for row in records)


def test_existing_rank_exclusions_remain_visible():
    records, summary = binary_score_subject_influence(example_components([1., 2., 3., 4., 0.]))
    assert summary["n_informative_subjects"] == 4
    assert not records[-1]["retained"] and records[-1]["precision_share"] == 0.
    assert np.isnan(records[-1]["pseudo_effect"])


@pytest.mark.parametrize("scale", [1e-6, 1e6])
def test_diagnostic_shares_are_invariant_to_target_units(scale):
    original = example_components([1., 2., 3., 4., 5.])
    transformed = example_components([1., 2., 3., 4., 5.])
    transformed.scores /= scale
    transformed.information /= scale**2
    transformed.reference_information /= scale**2
    transformed.biological_shapes *= scale**2
    records, summary = binary_score_subject_influence(original)
    transformed_records, transformed_summary = binary_score_subject_influence(transformed)
    np.testing.assert_allclose([row["precision_share"] for row in transformed_records], [row["precision_share"] for row in records], rtol=1e-7)
    assert transformed_summary["effective_weighted_subjects"] == pytest.approx(summary["effective_weighted_subjects"], rel=1e-7)
    assert transformed_summary["fitted_mean"] == pytest.approx(summary["fitted_mean"] * scale, rel=1e-7)


def test_duplicate_subjects_and_nonbinary_components_are_rejected():
    components = example_components([1.] * 6)
    components.subject_ids[1] = components.subject_ids[0]
    with pytest.raises(ValueError, match="unique aligned binary"):
        binary_score_subject_influence(components)
    components = example_components([1.] * 6)
    components.scores = np.ones((6, 2))
    with pytest.raises(ValueError, match="binary"):
        binary_score_subject_influence(components)


@pytest.mark.parametrize("inclusion", [.5, .03, 1e-7, 1 - 1e-7])
def test_archived_coordinate_transform_matches_direct_ec_score(inclusion):
    p = inclusion
    theta = np.array([[p * .2, p * .4, (1-p) * .6, .4], [p * .6, p * .2, (1-p) * .8, .2]])
    mapping = np.vstack([np.eye(4), [1, 1, 0, 0], [1, 0, 1, 1]])
    second = mapping * np.array([1., 2., 3., 1.])[None, :]
    counts = (np.array([[40., 17., 55., 60., 19., 4.], [15., 29., 13., 41., 9., 10.]]), np.array([[22., 24., 15., 20., 11., 14.], [33., 14., 19., 34., 18., 20.]]))
    original = efficient_shared_path_score(counts, (mapping, second), theta, [0, 0, 1, -1], [0, 1], 2, score_coordinate="ilr", return_reference=True)
    direct = efficient_shared_path_score(counts, (mapping, second), theta, [0, 0, 1, -1], [0, 1], 2, score_coordinate="proportion", return_reference=True)
    components = PathScoreComponents(original[0][None, :], original[1][None, :, :], original[2][None, :, :], np.array(["subject"]), (0, 1), [], [], "ilr", original[3][None, :, :])
    transformed = binary_score_components_to_proportions(components)
    for actual, expected in zip((transformed.scores[0], transformed.information[0], transformed.biological_shapes[0], transformed.reference_information[0]), direct):
        np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-10)
    np.testing.assert_array_equal(components.scores[0], original[0])
    assert components.score_coordinate == "ilr" and transformed.score_coordinate == "proportion"


def test_coordinate_transform_rejects_other_geometry_and_repeat_conversion():
    components = example_components([1.] * 5)
    with pytest.raises(ValueError, match="at least four"):
        binary_score_components_to_proportions(components)
    components.biological_shapes *= 4
    converted = binary_score_components_to_proportions(components)
    with pytest.raises(ValueError, match="binary ILR"):
        binary_score_components_to_proportions(converted)
