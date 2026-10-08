import numpy as np
import pytest
from tealeaf.sc.replication_audit import aligned_direction, coverage_correlation, reexpress_event_directions


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


def test_junction_directions_use_only_manifest_samples(tmp_path, monkeypatch):
    import json
    from types import SimpleNamespace
    import pandas as pd
    from scipy.sparse import csr_matrix
    import extra_scripts.audit_split_coverage_direction as audit

    folder = tmp_path / "reproducibility/fold0"
    folder.mkdir(parents=True)
    (folder / "contrasts.json").write_text(json.dumps([{"contrast_id": "contrast", "samples_a": ["a0"], "samples_b": ["b0"]}]))
    table = pd.DataFrame({"method": ["scQuint", "scQuint"], "contrast_id": ["contrast"] * 2, "feature_id": ["event", "missing"]})
    bundle = SimpleNamespace(counts=csr_matrix([[9, 1], [1, 9], [0, 10], [10, 0]]), samples=pd.DataFrame({"sample_id": ["a0", "b0", "a1", "b1"]}))
    monkeypatch.setattr(audit, "scquint_groups", lambda *args: {("scQuint", "contrast", "event"): {"contrast_id": "contrast", "indices": np.array([0, 1])}})
    result = audit.signed_junctions(table, bundle, tmp_path, 0, "scQuint")
    assert np.allclose(result.effect_vector.iloc[0], [-.8, .8])
    assert result.effect_features.iloc[0] == ["0", "1"]
    assert result.effect_vector.iloc[1] == []
    assert len(result) == len(table)


def test_majiq_directions_preserve_edge_identity(tmp_path):
    import pandas as pd
    from extra_scripts.audit_split_coverage_direction import signed_junctions
    from extra_scripts.assess_tilgner_junction_replication import majiq_feature_ids

    raw = pd.DataFrame({"gene_id": ["g", "g"], "seqid": ["chr1"] * 2, "start": [10, 10], "end": [20, 30], "a-raw_psi_quantile_0.500": [.2, .8], "b-raw_psi_quantile_0.500": [.6, .4]})
    folder = tmp_path / "reproducibility/fold0/majiq_min3_cov3/tests/raw"
    folder.mkdir(parents=True)
    raw.to_csv(folder / "contrast.tsv", sep="\t", index=False)
    features = majiq_feature_ids(raw)
    table = pd.DataFrame({"contrast_id": ["contrast"] * 3, "level_a": ["a"] * 3, "level_b": ["b"] * 3, "feature_id": [features.iloc[1], features.iloc[0], "missing"]})
    result = signed_junctions(table, None, tmp_path, 0, "MAJIQ Heterogen")
    assert np.allclose(result.effect_vector.iloc[:2].tolist(), [[-.4], [.4]])
    assert np.isnan(result.effect_vector.iloc[2][0])
    assert result.effect_features.iloc[0] == [features.iloc[1]]


def test_pair_completeness_never_averages_away_failed_subject_fits():
    import pandas as pd
    from tealeaf.sc.replication_audit import complete_paired_fits

    table = pd.DataFrame({"n_samples": [10, 10, 10, 10, 9], "n_subjects": [5, 4, 5, 5, 4], "converged": [True, True, True, False, True], "report_n_subjects": [5, 4, 4, 5, 4]})
    tested, reported = complete_paired_fits(table)
    assert tested.tolist() == [True, False, True, False, False]
    assert reported.tolist() == [True, False, False, False, False]


def test_effect_direction_ablation_preserves_family_ranking_and_external_zeros():
    import pandas as pd

    mapped = pd.DataFrame(dict(feature_id=list("abcde"), contrast_id=["A__B"] * 5, short_read_effect=[.2, -.3, .1, .2, .4], long_read_effect=[.3, .4, 0, .5, .6], replicate_1_dot_product=[.04, -.06, 0, .1, .2], replicate_2_dot_product=[.06, .09, 0, .1, .2], p_value=[.01, .02, .03, .04, .05], minimum_pooled_depth=[50] * 5))
    effects = pd.DataFrame(dict(feature_id=list("edcba"), contrast_id=["A__B"] * 5, score=[np.nan, 0, 1., 2., -3.]))
    result = reexpress_event_directions(mapped, effects, "score", "score direction")
    assert result.feature_id.tolist() == mapped.feature_id.tolist()
    np.testing.assert_array_equal(result.p_value, mapped.p_value)
    np.testing.assert_array_equal(result.minimum_pooled_depth, mapped.minimum_pooled_depth)
    assert result.pooled_replicated.tolist() == [False, True, False, False, False]
    assert result.both_replicates_replicated.tolist() == [False, False, False, False, False]
    assert result.direction_available.tolist() == [True, True, True, False, False]
    assert len(result) == len(mapped)
    assert mapped.short_read_effect.iloc[0] == .2
    with pytest.raises(ValueError, match="every fixed mapped"):
        reexpress_event_directions(mapped, effects.iloc[:-1], "score", "test")
    with pytest.raises(ValueError, match="unique identities"):
        reexpress_event_directions(mapped, pd.concat([effects, effects.iloc[:1]]), "score", "test")
