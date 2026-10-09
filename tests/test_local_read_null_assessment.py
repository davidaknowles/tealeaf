import pandas as pd
import pytest

from extra_scripts.assess_local_read_null import validate_trials, integration_recipe


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


def test_null_collation_preserves_different_driver_manifest_node_fields():
    assert integration_recipe(dict(quadrature_nodes=21), 'conditional') == dict(nodes=21, adaptive_integration=False)
    assert integration_recipe(dict(nodes=11, adaptive_integration=True), 'unconditional') == dict(nodes=11, adaptive_integration=True)


@pytest.mark.parametrize('manifest', [dict(nodes=True), dict(nodes=2), dict(nodes=11, adaptive_integration='true')])
def test_null_collation_rejects_undeclared_integration_types(manifest):
    with pytest.raises(ValueError):
        integration_recipe(manifest, 'unconditional')
