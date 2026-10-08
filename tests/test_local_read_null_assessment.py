import pandas as pd
import pytest

from extra_scripts.assess_local_read_null import validate_trials


def trials():
    return pd.DataFrame([dict(scenario=0, law=law, draw=draw, p_value=.5 if draw else 1., F_reference_p_value=.6 if draw else 1., converged=bool(draw)) for law in ("unconditional binomial", "conditional working model") for draw in range(2)])


def test_complete_trials_retain_unavailable_draws():
    validate_trials(trials(), 0, 2)


@pytest.mark.parametrize("problem", ["missing", "duplicate", "failed_p", "nonfinite", "wrong_scenario"])
def test_changed_or_incomplete_null_panels_are_rejected(problem):
    table = trials()
    if problem == "missing":
        table = table.iloc[:-1]
    elif problem == "duplicate":
        table = pd.concat([table, table.iloc[:1]])
    elif problem == "failed_p":
        table.loc[0, "p_value"] = .01
    elif problem == "nonfinite":
        table.loc[1, "p_value"] = float("nan")
    else:
        table.loc[0, "scenario"] = 1
    with pytest.raises(ValueError):
        validate_trials(table, 0, 2)
