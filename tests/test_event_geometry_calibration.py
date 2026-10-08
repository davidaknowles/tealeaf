import numpy as np
import pandas as pd
import pytest

from extra_scripts.assess_event_geometry_calibration import calibrate_tables, cohort, count_null, replace_fixed_mapping_probabilities
from test_event_score_subject_influence import archive_fixture
from types import SimpleNamespace
import json
from tealeaf.sc.path_score_mixed import MODEL_VERSION


def fixture():
    observed = pd.DataFrame(dict(test_id=list("abcdx"), n_subjects=[5] * 4 + [0], raw_p_value=[.001, .01, .001, .01, 1.], p_value=[.1] * 4 + [1.], calibration_stratum=["1|5"] * 4 + ["1|0"], complete_subject_fits=[True] * 4 + [False]))
    geometry = pd.DataFrame(dict(test_id=list("abcd"), n_informative_subjects=[5] * 4, geometry_class=["balanced"] * 2 + ["dominated"] * 2))
    null = pd.DataFrame([dict(test_id=key, replicate=draw, raw_p_value=.2 if key in "ab" else .0001, p_value=.5) for key in "abcd" for draw in range(32)])
    return observed, null, geometry


def test_calibration_never_mixes_geometry_or_removes_failed_tests():
    observed, null, geometry = fixture()
    before, before_null = observed.copy(deep=True), null.copy(deep=True)
    result, draws = calibrate_tables(observed, null, geometry)
    np.testing.assert_array_equal(result.test_id, observed.test_id)
    assert result.loc[result.test_id.eq("x"), "p_value"].item() == 1.
    assert result.loc[result.test_id.isin(list("ab")), "p_value"].eq(1 / 33).all()
    assert result.loc[result.test_id.isin(list("cd")), "p_value"].eq(1.).all()
    assert len(draws) == 128 and len(result) == 5
    np.testing.assert_array_equal(result.raw_p_value, observed.raw_p_value)
    pd.testing.assert_frame_equal(observed, before)
    pd.testing.assert_frame_equal(null, before_null)


def test_observed_and_geometry_row_order_cannot_change_labels():
    observed, null, geometry = fixture()
    original, _ = calibrate_tables(observed, null, geometry)
    shuffled, _ = calibrate_tables(observed.sample(frac=1., random_state=2), null.sample(frac=1., random_state=3), geometry.sample(frac=1., random_state=4))
    pd.testing.assert_frame_equal(original.set_index("test_id").sort_index(), shuffled.set_index("test_id").sort_index())


def prepare_cohort(tmp_path):
    root, merged = archive_fixture(tmp_path)
    shard = root / "shard_0"
    (shard / "failures.json").write_text("[]\n")
    pd.DataFrame(dict(test_id=list("ab"), n_expected_subjects=[5, 5], score_coordinate=["ilr", "ilr"])).to_csv(shard / "score_contexts.tsv.gz", sep="\t", index=False)
    table = pd.read_csv(merged / "paired_path.tsv", sep="\t")
    table["calibration_stratum"] = "1|5"
    table.to_csv(merged / "paired_path.tsv", sep="\t", index=False)
    null = pd.read_csv(shard / "paired_path_null.tsv.gz", sep="\t")
    null["raw_p_value"] = null.p_value
    null.to_csv(merged / "paired_path_null.tsv.gz", sep="\t", index=False)
    output = tmp_path / "geometry"
    args = SimpleNamespace(cohort_root=root, merged_dir=merged, output_dir=output, shard_count=1)
    cohort(args)
    return args


def test_complete_cohort_pipeline_and_declared_context_guard(tmp_path):
    args = prepare_cohort(tmp_path)
    output = args.output_dir
    result = pd.read_csv(output / "paired_path.tsv", sep="\t").set_index("test_id")
    assert result.at["a", "geometry_class"] == "dominated"
    assert result.at["b", "geometry_class"] == "balanced"
    assert len(result) == 2 and result.p_value.eq(1.).all()
    with pytest.raises(ValueError, match="new output"):
        cohort(args)


@pytest.mark.parametrize("failed,missing", [(False, False), (True, False), (False, True)])
def test_count_null_pipeline_preserves_failed_draws_and_requires_whole_family(tmp_path, failed, missing):
    prepared = prepare_cohort(tmp_path)
    source = prepared.cohort_root / "shard_0"
    null_root = tmp_path / "count_null"
    shard = null_root / "shard_0"
    shard.mkdir(parents=True)
    table = pd.read_csv(source / "paired_path.tsv", sep="\t")
    table["draw"], table["converged"] = 0, True
    table["n_expected_subjects"], table["n_fitted_subjects"] = 5, 5
    if failed:
        table.loc[table.test_id.eq("a"), "converged"] = False
    if missing:
        table = table.iloc[:1]
    table.to_csv(shard / "observed.tsv", sep="\t", index=False)
    diagnostics = pd.read_csv(source / "subject_scores.tsv.gz", sep="\t")
    diagnostics["draw"] = 0
    diagnostics.to_csv(shard / "subject_null_diagnostics.tsv.gz", sep="\t", index=False)
    settings = dict(requested_ids=list("ab"), draws=1, expected_strategies=["fixture"], mixed_score_version=MODEL_VERSION, score_coordinate="ilr", information_metric="reference")
    (shard / "settings.json").write_text(json.dumps(settings))
    args = SimpleNamespace(training_dir=prepared.output_dir, null_root=null_root, shard_count=1, output_dir=tmp_path / "assessed_null")
    if missing:
        with pytest.raises(ValueError, match="missing, duplicated or foreign"):
            count_null(args)
        return
    count_null(args)
    result = pd.read_csv(args.output_dir / "trials.tsv.gz", sep="\t").set_index("test_id")
    assert len(result) == 2 and set(result.index) == set("ab")
    assert result.geometry_calibrated_p_value.eq(1.).all()
    if failed:
        assert result.at["a", "geometry_class"] == "failed"
        assert result.at["a", "p_value"] == 1.
    receipt = json.loads((args.output_dir / "manifest.json").read_text())
    assert receipt["complete_subject_fits"] == (1 if failed else 2)
    assert "own real parent excluded" in receipt["calibration"]


def test_fixed_long_read_mapping_changes_only_calibration_not_selected_associations():
    mapped = pd.DataFrame(dict(feature_id=["a", "b"], contrast_id=["x", "y"], p_value=[.01, .2], raw_p_value=[.0001, .1], long_read_effect=[-.1, 0.], short_read_effect=[.2, .3], fdr=[.03, .3]))
    tests = pd.DataFrame(dict(feature_id=["b", "a", "extra"], contrast_id=["y", "x", "z"], p_value=[.9, .0003, .1], raw_p_value=[.1, .0001, .01], legacy_calibrated_p_value=[.2, .01, .2], fdr=[1., .0005, .2]))
    before = mapped.copy(deep=True)
    result = replace_fixed_mapping_probabilities(mapped, tests)
    pd.testing.assert_frame_equal(result.drop(columns=["p_value", "fdr"]), mapped.drop(columns=["p_value", "fdr"]))
    np.testing.assert_array_equal(result.p_value, [.0003, .9])
    pd.testing.assert_frame_equal(mapped, before)
    for changed in (tests.loc[tests.feature_id.ne("a")], tests.assign(raw_p_value=.5)):
        with pytest.raises(ValueError, match="missing tests or differs"):
            replace_fixed_mapping_probabilities(mapped, changed)


@pytest.mark.parametrize("problem", ["missing_geometry", "unknown_class", "missing_draw", "duplicate_draw", "wrong_rank"])
def test_incomplete_geometry_or_original_null_is_rejected(problem):
    observed, null, geometry = fixture()
    if problem == "missing_geometry":
        geometry = geometry.iloc[:-1]
    elif problem == "unknown_class":
        geometry.loc[0, "geometry_class"] = "favorable"
    elif problem == "missing_draw":
        null = null.iloc[:-1]
    elif problem == "duplicate_draw":
        null = pd.concat([null, null.iloc[:1]], ignore_index=True)
    else:
        geometry.loc[0, "n_informative_subjects"] = 4
    with pytest.raises(ValueError):
        calibrate_tables(observed, null, geometry)
