import pandas as pd
import pytest

from extra_scripts.compare_count_null_implementations import compare_tables, normalize_settings


def test_numerical_comparison_preserves_failed_trials_and_tail_decisions():
    table = pd.DataFrame(dict(test_id=["a", "b"], draw=[0, 0], p_value=[.001, 1.], statistic=[10., 0.]))
    changed = table.iloc[::-1].copy()
    changed.loc[changed.test_id.eq("a"), "p_value"] += 1e-10
    result = compare_tables(table, changed, ["test_id", "draw"], ["p_value", "statistic"])
    assert result["trials"] == 2 and result["p_value_matches"]
    # A tolerance agreement does not excuse crossing a diagnostic threshold.
    assert result["p_decisions_differ_at_0.001"] == 0
    changed.loc[changed.test_id.eq("a"), "p_value"] -= 2e-10
    assert compare_tables(table, changed, ["test_id", "draw"], ["p_value"])["p_decisions_differ_at_0.001"] == 1
    with pytest.raises(ValueError, match="trial family"):
        compare_tables(table, changed.iloc[:1], ["test_id", "draw"], ["p_value"])


def test_known_legacy_defaults_do_not_discard_opportunity_or_coordinate_changes():
    old = dict(information_metric="reference", seed=15)
    explicit = dict(old, count_likelihood="multinomial", ec_opportunity_scale=0., kernel_units="prepared", simulation_kernel_units="analysis", scalar_fast=False)
    assert normalize_settings(old) == normalize_settings(explicit)
    assert normalize_settings(old) != normalize_settings(dict(explicit, kernel_units="fragment"))
    assert normalize_settings(old) != normalize_settings(dict(explicit, ec_opportunity_scale=1.))
