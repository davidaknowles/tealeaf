import numpy as np
import pandas as pd
import pytest

from extra_scripts.audit_lr_fit_availability import audit_prefix


def fixture():
    prefix = pd.DataFrame(dict(rank=range(1, 201), feature_id=[f"e{i}" for i in range(200)], contrast_id="cell_type__A__B"))
    fits = prefix.iloc[:199].drop(columns="rank").assign(converged=True, complete_subject_fits=True, complete_reporting_fits=True, n_subjects=8, n_expected_subjects=8, n_fitted_subjects=8, p_value=.01, raw_p_value=.01, effect_size=.2, test_ilr_effect_size=.4)
    return prefix, fits


def test_incomplete_and_missing_fits_remain_in_fixed_rank_prefix():
    prefix, fits = fixture()
    fits.loc[1, ["converged", "complete_subject_fits", "complete_reporting_fits"]] = False
    fits.loc[2, "complete_reporting_fits"] = False
    fits.loc[3, "effect_size"] = np.nan
    fits.loc[4, "test_ilr_effect_size"] = 0.
    table = audit_prefix(prefix, fits)
    assert len(table) == 200
    assert table.requested_in_completed_candidate.sum() == 199
    assert table.complete_subject_fits.sum() == 198
    assert table.finite_usage_direction.sum() == 196
    assert table.finite_score_direction.sum() == 197


def test_duplicate_or_incomplete_native_prefix_is_not_a_valid_audit():
    prefix, fits = fixture()
    with pytest.raises(ValueError, match="top-200"):
        audit_prefix(prefix.iloc[:100], fits)
    with pytest.raises(ValueError, match="unique"):
        audit_prefix(prefix, pd.concat([fits, fits.iloc[:1]]))
