import pandas as pd

from extra_scripts.plot_tilgner_method_replication import rank_agreement_table


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
    assert observed.groupby("method").size().to_dict() == {"MAJIQ Heterogen": 2, "Tealeaf": 2}
    final = observed.groupby("method")["cumulative_agreement"].last().to_dict()
    assert final == {"MAJIQ Heterogen": 0.5, "Tealeaf": 0.5}
