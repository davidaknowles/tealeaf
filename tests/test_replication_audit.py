import numpy as np
import pytest
from tealeaf.sc.replication_audit import aligned_direction, coverage_correlation


def test_coverage_uses_p_not_negative_log_p_and_filters_missing():
    result = coverage_correlation([.8, .4, .1, np.nan], [1, 2, 3, 4])
    assert np.isclose(result["rho_p_coverage"], -1)
    assert result["n"] == 3
    assert np.isnan(coverage_correlation([1, 1, 1], [1, 2, 3])["rho_p_coverage"])


def test_partial_correlation_detects_depth_fully_explained_by_control():
    result = coverage_correlation(np.arange(10)[::-1] / 10, np.arange(1, 11), np.arange(10))
    assert np.isnan(result["partial_rho"])


def test_direction_alignment_is_feature_not_column_based():
    result = aligned_direction([1, -1], [-2, 2], ["a", "b"], ["b", "a"])
    assert result["direction_agrees"]
    assert result["agreeing_components"] == 2
    assert np.isclose(result["cosine"], 1)
    assert not aligned_direction([1], [-1])["direction_agrees"]
    assert np.isnan(aligned_direction([0], [1])["direction_agrees"])
    with pytest.raises(ValueError, match="sets differ"):
        aligned_direction([1], [1], ["a"], ["b"])


def test_gene_coverage_does_not_weight_pairs_by_number_of_blocks():
    import pandas as pd
    from extra_scripts.audit_split_coverage_direction import gene_table
    table = pd.DataFrame({"gene_id": ["g"] * 11, "pair_id": ["a||b"] * 10 + ["b||c"], "p_value": [.1] * 11, "coverage": [10.] * 10 + [30.], "reference_n_subjects": [8.] * 11})
    result = gene_table(table)
    assert result.coverage.iloc[0] == 20
    assert result.n_features.iloc[0] == 11


def test_direction_union_is_symmetric_and_uses_consistent_level_order():
    import pandas as pd
    from extra_scripts.audit_split_coverage_direction import direction_rows, gene_table
    first = pd.DataFrame({"gene_id": ["g", "h"], "pair_id": ["a||b"] * 2, "feature_id": ["e1", "e2"], "p_value": [.001, .5], "coverage": [10, 20], "reference_n_subjects": [8, 8], "level_a": ["a", "a"], "level_b": ["b", "b"], "effect_size": [1., -1.]})
    second = first.copy()
    second["p_value"] = [.5, .001]
    second["level_a"], second["level_b"] = "b", "a"
    second["effect_size"] = [-1., 1.]
    tables = [first, second]
    records, _ = direction_rows(tables, [gene_table(table) for table in tables], "test", "test")
    assert all(record["event_BH_union"] for record in records)
    assert not any(record["event_BH_intersection"] for record in records)
    assert all(record["direction_agrees"] for record in records)
