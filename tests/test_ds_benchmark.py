import numpy as np
import pandas as pd

from tealeaf.sc.ds_benchmark import (
    aggregate_feature_pvalues,
    aggregate_gene_pvalues,
    shared_gene_reproducibility,
    simes_pvalue,
)


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
