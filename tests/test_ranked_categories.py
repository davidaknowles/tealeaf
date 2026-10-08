import numpy as np
import pandas as pd
import pytest

from tealeaf.sc.replication_audit import ranked_category_summary


def fixture():
    return pd.DataFrame(dict(method=["a"] * 5, rank=[1, 2, 3, 4, 5], event_type=["AF", "AF", "SE", None, "SE"], pooled_replicated=[False, True, True, False, True]))


def test_categories_use_global_prefix_without_reranking_or_dropping_zeros():
    rows = pd.DataFrame(ranked_category_summary(fixture(), cutoffs=(3, 5)))
    top3 = rows.loc[rows.cutoff.eq(3)].set_index("category")
    assert set(top3.index) == {"AF", "SE"}
    assert top3.loc["AF", "n_tests"] == 2
    assert top3.loc["AF", "agreement"] == .5
    assert top3.loc["SE", "n_tests"] == 1
    assert np.isclose(top3.fraction_of_prefix.sum(), 1)
    assert rows.loc[rows.category.eq("unclassified"), "n_tests"].iloc[0] == 1


def test_short_prefix_is_explicit_not_extrapolated():
    rows = pd.DataFrame(ranked_category_summary(fixture(), cutoffs=(10,)))
    assert not rows.complete_prefix.any()
    assert rows.n_tests.sum() == 5


@pytest.mark.parametrize("defect", ["missing_category", "rank_gap", "invalid_boolean"])
def test_invalid_rank_category_inputs_fail(defect):
    table = fixture()
    if defect == "missing_category":
        table = table.drop(columns="event_type")
    elif defect == "rank_gap":
        table.loc[0, "rank"] = 6
    else:
        table["pooled_replicated"] = table.pooled_replicated.astype(object)
        table.loc[0, "pooled_replicated"] = "unknown"
    with pytest.raises(ValueError):
        ranked_category_summary(table)
