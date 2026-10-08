import numpy as np
import pandas as pd
import pytest

from extra_scripts.summarize_ec_count_null import validate_event_mass_truth, validate_requested_trials


def test_requested_null_trials_retain_failures_and_reject_missing_keys():
    settings = {"requested_ids": ["a", "b"], "draws": 2, "expected_strategies": ["x", "y"]}
    table = pd.DataFrame([{"test_id": test_id, "draw": draw, "strategy": strategy, "p_value": .001 if test_id == "a" else np.nan, "converged": test_id == "a"} for test_id in ("a", "b") for draw in range(2) for strategy in ("x", "y")])
    result = validate_requested_trials(table, settings)
    assert len(result) == 8
    assert result.loc[result.test_id.eq("b"), "p_value"].eq(1).all()
    assert result.loc[result.test_id.eq("a"), "p_value"].eq(.001).all()
    with pytest.raises(ValueError, match="missing or duplicate"):
        validate_requested_trials(table.iloc[:-1], settings)
    with pytest.raises(ValueError, match="missing or duplicate"):
        validate_requested_trials(pd.concat([table, table.iloc[:1]]), settings)
    table.loc[0, "p_value"] = np.nan
    with pytest.raises(ValueError, match="invalid successful"):
        validate_requested_trials(table, settings)


def mass_fixture():
    settings = dict(requested_ids=["a", "b", "failed"], draws=1, residual_concentration=None)
    observed = pd.DataFrame(dict(test_id=["a", "b", "failed"], draw=[0, 0, 0], converged=[True, True, False]))
    truth = pd.DataFrame(dict(test_id=["a", "b"], draw=[0, 0], has_outside_transcripts=["True", "False"], tilt_applied=["True", "False"], maximum_absolute_subject_conditional_path_change=[0., 1e-16]))
    return truth, observed, settings


def test_mass_truth_distinguishes_generated_and_failed_trials_without_losing_denominators():
    truth, observed, settings = mass_fixture()
    result = validate_event_mass_truth(truth, observed, settings)
    assert result.tilt_applied.tolist() == [True, False]
    assert len(observed) == 3 and len(result) == 2
    settings["residual_concentration"] = 20
    truth.loc[0, "maximum_absolute_subject_conditional_path_change"] = .3
    assert len(validate_event_mass_truth(truth, observed, settings)) == 2


@pytest.mark.parametrize("defect", ["missing", "duplicate", "foreign", "wrong_tilt", "invalid_flag", "nonnull", "nonfinite"])
def test_invalid_mass_truth_cannot_be_used_as_target_null_evidence(defect):
    truth, observed, settings = mass_fixture()
    if defect == "missing":
        truth = truth.iloc[1:]
    elif defect == "duplicate":
        truth = pd.concat([truth, truth.iloc[:1]])
    elif defect == "foreign":
        truth.loc[0, "test_id"] = "foreign"
    elif defect == "wrong_tilt":
        truth.loc[0, "tilt_applied"] = "False"
    elif defect == "invalid_flag":
        truth.loc[0, "has_outside_transcripts"] = "unknown"
    elif defect == "nonnull":
        truth.loc[0, "maximum_absolute_subject_conditional_path_change"] = .01
    else:
        truth.loc[0, "maximum_absolute_subject_conditional_path_change"] = np.nan
    with pytest.raises(ValueError):
        validate_event_mass_truth(truth, observed, settings)
