import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from extra_scripts.audit_event_information import cross_ranking_effects, event_evidence, gene_support, rank_events, transcript_evidence


def mapping(method, effects, pvalues, truth):
    return pd.DataFrame({"method": method, "contrast_id": "A_B", "feature_id": [f"e{x}" for x in range(len(effects))], "p_value": pvalues, "raw_p_value": pvalues, "short_read_effect": effects, "long_read_effect": truth})


def test_ranking_direction_crossing_does_not_change_the_shared_event_family():
    native = mapping("native", np.ones(200), np.arange(1, 201) / 201, np.ones(200))
    hybrid = mapping("hybrid", np.r_[np.ones(100), -np.ones(100)], np.arange(200, 0, -1) / 201, np.ones(200))
    summary, shared, _, _ = cross_ranking_effects(native, hybrid)
    common = summary.loc[summary.scope.eq("shared event-contrast universe") & summary.cutoff.eq(100)]
    assert len(shared) == 200
    assert len(common) == 4
    assert set(common.n_available) == {200}
    assert common.loc[common.ranking.eq("native_p") & common.direction.eq("hybrid_effect"), "agreement"].iloc[0] == 1
    assert common.loc[common.ranking.eq("hybrid_p") & common.direction.eq("hybrid_effect"), "agreement"].iloc[0] == 0


def test_crossing_rejects_external_truth_or_duplicate_identity_mismatch():
    first = mapping("native", [1], [.01], [.3])
    second = mapping("hybrid", [1], [.01], [.2])
    with pytest.raises(ValueError, match="same LR"):
        cross_ranking_effects(first, second)
    with pytest.raises(ValueError, match="unique"):
        cross_ranking_effects(pd.concat([first, first]), first)


def test_gene_support_preserves_supported_zero_expression_transcripts():
    # No observed expression or external agreement appears in this decision.
    support = gene_support(["g.1"], [np.array([0, 1, 2])], [np.array([0, 1])], [sparse.csr_matrix([[1, 0, 0], [0, 1, 0]])], ["t1.1", "t2.2", "t3.1"])
    assert support["g"]["supported_transcripts"] == {"t1", "t2"}
    assert support["g"]["all_transcripts"] == {"t1", "t2", "t3"}


def test_zero_external_effect_policy_is_explicit_and_never_silently_dropped():
    first = mapping("native", [1, 1], [.01, .02], [0, .3])
    second = mapping("hybrid", [1, 1], [.01, .02], [0, .3])
    summary, shared, native, _ = cross_ranking_effects(first, second)
    assert len(shared) == len(native) == 2
    _, restricted, _, _ = cross_ranking_effects(first, second, exclude_zero_lr=True)
    assert len(restricted) == 1


def test_rank_ties_follow_published_raw_statistic_then_feature_order():
    table = pd.DataFrame({"feature_id": ["b", "a", "c"], "native_p": [.1, .1, .1], "native_p_raw": [.02, .02, .03], "native_p_statistic": [2, 1, 100], "pooled_replicated": ["True", "False", None]})
    assert rank_events(table, "native_p").feature_id.tolist() == ["b", "a", "c"]


def test_event_support_distinguishes_global_read_support_and_expression():
    represented, supported, mass = transcript_evidence([sparse.csr_matrix([[1, 0, 1]])], ["a.1", "b.2", "c.1"], sparse.csr_matrix([[2, 3, 5, 7]]), ["a.1", "b.1", "c.1", "d.1"])
    evidence = event_evidence({"a", "b"}, {"c", "d"}, {"a"}, represented, supported, mass)
    assert evidence["n_missing"] == 3
    assert evidence["n_missing_global_supported"] == 1
    assert evidence["n_missing_global_zero_support"] == 1
    assert evidence["n_missing_not_in_features"] == 1
    assert evidence["n_missing_positive_expression"] == 3
    assert evidence["missing_expression_fraction"] == 15 / 17
    assert evidence["included_missing_expression_fraction"] == 3 / 5
    assert evidence["excluded_missing_expression_fraction"] == 1
