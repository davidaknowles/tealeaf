import pandas as pd

from extra_scripts.plot_tilgner_method_replication import rank_agreement_table, replace_rank_methods
from tealeaf.sc.replication_audit import ranked_direction_summary


def test_rank_agreement_uses_only_eligible_directional_calls(tmp_path):
    comparator = pd.DataFrame(
        [
            ["MAJIQ Heterogen", "m1", True, True, 20, 0.01],
            ["MAJIQ Heterogen", "m2", True, False, 20, 0.02],
            ["MAJIQ Heterogen", "m3", True, None, 20, 0.03],
        ],
        columns=["method", "feature_id", "mapping_complete", "pooled_replicated", "minimum_pooled_depth", "p_value"],
    )
    tealeaf = pd.DataFrame(
        [["t1", True, True, 20, 0.01], ["t2", True, False, 20, 0.02]],
        columns=["test_id", "mapping_complete", "pooled_replicated", "minimum_pooled_depth", "p_value"],
    )
    comparator_path = tmp_path / "junction.tsv"
    tealeaf_path = tmp_path / "tealeaf.tsv"
    comparator.to_csv(comparator_path, sep="\t", index=False)
    tealeaf.to_csv(tealeaf_path, sep="\t", index=False)
    observed = rank_agreement_table(comparator_path, tealeaf_path, None)
    assert observed.groupby("method").size().to_dict() == {"MAJIQ Heterogen": 2, "Tealeaf pairwise": 2}
    final = observed.groupby("method")["cumulative_agreement"].last().to_dict()
    assert final == {"MAJIQ Heterogen": 0.5, "Tealeaf pairwise": 0.5}


def test_rank_agreement_breaks_calibrated_pvalue_ties(tmp_path):
    comparator = pd.DataFrame(
        [["MAJIQ Heterogen", f"m{i}", True, bool(i % 2), 20, 0.01]
         for i in range(3)],
        columns=["method", "feature_id", "mapping_complete", "pooled_replicated", "minimum_pooled_depth", "p_value"],
    )
    tealeaf = pd.DataFrame(
        [[f"t{i}", True, bool(i != 1), 20, 0.01] for i in range(3)],
        columns=["test_id", "mapping_complete", "pooled_replicated", "minimum_pooled_depth", "p_value"],
    )
    comparator_path = tmp_path / "junction.tsv"
    tealeaf_path = tmp_path / "tealeaf.tsv"
    comparator.to_csv(comparator_path, sep="\t", index=False)
    tealeaf.to_csv(tealeaf_path, sep="\t", index=False)
    observed = rank_agreement_table(comparator_path, tealeaf_path, None)
    tealeaf_rows = observed[observed.method.eq("Tealeaf pairwise")]
    assert tealeaf_rows.p_tie_size.eq(3).all()
    assert tealeaf_rows["rank"].tolist() == [1, 2, 3]


def test_all_tested_rank_replacement_keeps_nonsignificant_tests():
    import numpy as np
    replacement = pd.DataFrame({"method": ["LeafCutter"] * 210, "contrast_id": ["a_b"] * 210, "feature_id": [f"e{i}" for i in range(210)], "p_value": np.linspace(.1, .9, 210), "q_value": [1.] * 210, "mapping_complete": [True] * 210, "minimum_pooled_depth": [20] * 210, "pooled_replicated": [True] * 210})
    previous = replacement.iloc[:1].assign(rank=1)
    observed = replace_rank_methods(previous, replacement)
    assert len(observed) == 200
    assert observed["rank"].tolist() == list(range(1, 201))
    assert observed.q_value.eq(1).all()


def test_rank_area_weights_early_hits_and_does_not_extrapolate():
    import numpy as np
    table = pd.DataFrame({"method": ["a"] * 3, "rank": [1, 2, 3], "pooled_replicated": [True, False, True]})
    summary = ranked_direction_summary(table, (3, 100))
    assert np.isclose(summary[0]["agreement"], 2 / 3)
    assert np.isclose(summary[0]["normalized_auc"], (1 + .5 + 2 / 3) / 3)
    assert np.isnan(summary[1]["normalized_auc"])


def test_omnibus_ranking_uses_its_own_continuous_tail_not_pairwise_tail(tmp_path):
    empty = pd.DataFrame(columns=["method", "feature_id", "mapping_complete", "pooled_replicated", "minimum_pooled_depth", "p_value"])
    pairs = pd.DataFrame({"test_id": ["p1", "p2"], "block_id": ["b1", "b2"], "p_value": [.01, .01], "raw_p_value": [1e-12, .1], "statistic": [100., 1.], "original_effect_norm": [.2, .3], "mapping_complete": [True, True], "pooled_replicated": [False, True], "minimum_pooled_depth": [20, 20]})
    omnibus = pd.DataFrame({"block_id": ["b1", "b2"], "p_value": [.001, .001], "raw_p_value": [.1, 1e-12], "statistic": [1., 100.], "fdr": [.01, .01]})
    paths = [tmp_path / name for name in ("junction.tsv", "pairs.tsv", "omnibus.tsv")]
    for table, path in zip((empty, pairs, omnibus), paths):
        table.to_csv(path, sep="\t", index=False)
    ranked = rank_agreement_table(paths[0], paths[1], None, omnibus_path=paths[2])
    selected = ranked.loc[ranked.method.eq("Tealeaf omnibus")]
    assert selected.feature_id.tolist() == ["b2", "b1"]
    assert selected.pooled_replicated.tolist() == [True, False]
