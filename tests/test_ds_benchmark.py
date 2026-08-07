import numpy as np
import pandas as pd

from tealeaf.sc.ds_benchmark import (
    aggregate_feature_pvalues,
    aggregate_gene_pvalues,
    aggregate_gene_pair_pvalues,
    calibrate_pvalues_from_null,
    compositional_pairwise_table,
    filter_grouped_subject_records,
    shared_pair_gene_reproducibility,
    shared_gene_reproducibility,
    simes_pvalue,
)


def test_filter_grouped_subject_records_drops_other_subjects():
    grouped = {
        "b1": [{"mouse": "m1", "x": 1}, {"mouse": "m2", "x": 2}],
        "b2": [{"mouse": "m3", "x": 3}],
    }
    result = filter_grouped_subject_records(grouped, {"m2"})
    assert result == {"b1": [{"mouse": "m2", "x": 2}]}


def test_compositional_pairwise_table_normalizes_and_filters():
    table = pd.DataFrame({
        "gene_id": ["g1.2", "g2"],
        "contrast": ["a_vs_b", "a_vs_c"],
        "p_value": [0.01, 0.02],
        "converged": [True, True],
        "inference_eligible": [True, False],
        "n_mice": [8, 9],
    })
    result = compositional_pairwise_table(
        table, method="Tealeaf DM2", min_paired_subjects=8
    )
    assert result.to_dict("records") == [{
        "method": "Tealeaf DM2",
        "gene_id": "g1.2",
        "p_value": 0.01,
        "level_a": "a",
        "level_b": "b",
    }]


def test_calibrate_pvalues_from_null_uses_method_specific_ecdf():
    table = pd.DataFrame({
        "method": ["a", "a", "b"],
        "p_value": [0.01, 0.5, 0.01],
    })
    null = pd.DataFrame({
        "method": ["a", "a", "a", "b"],
        "p_value": [0.001, 0.02, 0.8, 0.5],
    })
    result = calibrate_pvalues_from_null(table, null)
    np.testing.assert_allclose(result.p_value, [0.5, 0.75, 0.5])
    np.testing.assert_allclose(result.raw_p_value, table.p_value)


def test_simes_pvalue():
    assert simes_pvalue([0.01, 0.04, 0.5]) == 0.03


def test_aggregate_gene_pvalues_strips_versions():
    table = pd.DataFrame(
        {
            "method": ["a", "a", "a"],
            "gene_id": ["g1.2", "g1.2", "g2"],
            "p_value": [0.01, 0.2, 0.5],
        }
    )
    result = aggregate_gene_pvalues(table).set_index("gene_id")
    assert result.loc["g1", "p_value"] == 0.02
    assert result.loc["g1", "n_features"] == 2


def test_aggregate_feature_pvalues_collapses_pairwise_tests_first():
    table = pd.DataFrame({
        "method": ["m"] * 4,
        "feature_id": ["b1", "b1", "b1", "b2"],
        "gene_id": ["g1.2"] * 4,
        "p_value": [0.01, 0.2, 0.3, 0.04],
    })
    result = aggregate_feature_pvalues(table).set_index("feature_id")
    assert np.isclose(result.loc["b1", "p_value"], 0.03)
    assert result.loc["b1", "n_pairwise"] == 3
    assert result.loc["b1", "gene_id"] == "g1"

    gene = aggregate_gene_pvalues(result.reset_index())
    assert np.isclose(gene.loc[0, "p_value"], 0.04)


def test_shared_gene_reproducibility_uses_cross_method_intersection():
    fold0 = pd.DataFrame(
        {
            "method": ["a", "a", "b", "b"],
            "gene_id": ["g1", "g2", "g1", "g3"],
            "p_value": [1e-8, 0.5, 1e-7, 0.4],
            "n_features": 1,
        }
    )
    fold1 = pd.DataFrame(
        {
            "method": ["a", "a", "b", "b"],
            "gene_id": ["g1", "g4", "g1", "g5"],
            "p_value": [1e-9, 0.5, 1e-6, 0.4],
            "n_features": 1,
        }
    )
    metrics, topk, genes = shared_gene_reproducibility(
        [fold0, fold1], top_k=(1,), reference_method="a"
    )
    assert metrics.shared_genes.tolist() == [1, 1]
    assert metrics.replicated_bh.tolist() == [1, 1]
    assert topk.overlap_fraction.tolist() == [1.0, 1.0]
    assert set(genes.gene_id) == {"g1"}


def test_pairwise_reproducibility_matches_gene_and_level_pair():
    fold0 = pd.DataFrame({
        "method": ["a", "a", "b", "b"],
        "gene_id": ["g1", "g1", "g1", "g1"],
        "level_a": ["x", "x", "x", "x"],
        "level_b": ["y", "z", "y", "w"],
        "p_value": [1e-8, 1e-10, 1e-7, 1e-12],
    })
    fold1 = fold0.copy()
    paired = aggregate_gene_pair_pvalues(fold0)
    assert set(paired.pair_id) == {"x||y", "x||z", "w||x"}
    metrics, _, genes = shared_pair_gene_reproducibility(
        [fold0, fold1], top_k=(1,), reference_method="a"
    )
    assert metrics.shared_gene_pairs.tolist() == [1, 1]
    assert metrics.replicated_bh.tolist() == [1, 1]
    assert set(genes.gene_id) == {"g1"}
