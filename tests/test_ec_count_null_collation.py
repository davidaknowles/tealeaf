import numpy as np
import pandas as pd
import pytest

from extra_scripts.summarize_ec_count_null import validate_requested_trials


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
