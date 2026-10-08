import numpy as np
import pandas as pd
import pytest

from extra_scripts.audit_event_null_coverage import attach_frozen_coverage, coverage_summaries


def fixture():
    attributes = pd.DataFrame(dict(test_id=["low", "high"], median_gene_umis=[25., 10000.], expected_subjects=[4, 12], transcripts=[3, 10], ecs=[5, 30], coverage_quartile=[1, 4]))
    trials = pd.DataFrame(dict(test_id=["low", "low", "high", "high"], draw=[0, 1, 0, 1], strategy=["score"] * 4, scenario=["strict"] * 4, converged=["False", "True", "True", "True"], p_value=[1., .1, .02, .001], raw_p_value=[1., .2, .03, .002]))
    return trials, attributes


def test_failed_trials_keep_prefit_coverage_and_all_denominators():
    trials, attributes = fixture()
    table = attach_frozen_coverage(trials, attributes)
    assert len(table) == 4 and table.median_gene_umis.iloc[0] == 25
    assert not table.converged.iloc[0]
    strata, correlations = coverage_summaries(table)
    assert strata.requested_trials.tolist() == [2, 2]
    assert strata.fitted_trials.tolist() == [1, 2]
    np.testing.assert_array_equal(strata["calibrated_rate_0.05"], [0, 1])
    assert set(correlations.scope) == {"all requested trials", "successful-fit rho only, not calibration denominator"}
    assert set(correlations.draw) == {"all repeated draws", "0", "1"}


@pytest.mark.parametrize("defect", ["missing", "duplicate", "failed_p", "zero_coverage"])
def test_invalid_null_coverage_family_is_rejected(defect):
    trials, attributes = fixture()
    if defect == "missing":
        attributes = attributes.iloc[:1]
    elif defect == "duplicate":
        attributes = pd.concat([attributes, attributes.iloc[:1]])
    elif defect == "failed_p":
        trials.loc[0, "p_value"] = .001
    else:
        attributes.loc[0, "median_gene_umis"] = 0
    with pytest.raises(ValueError):
        attach_frozen_coverage(trials, attributes)
