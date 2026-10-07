import numpy as np
import pytest

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_pooling import quantify_effective_paths, joint_path_dm_test


def direct_data(subjects=12, paths=2):
    labels = np.tile([0, 1], subjects)
    clusters = np.repeat(np.arange(subjects), 2)
    rng = np.random.default_rng(91)
    proportions = rng.dirichlet(np.ones(paths) * 10, size=subjects)[:, None, :].repeat(2, axis=1).reshape(-1, paths)
    depths = np.tile([200., 2000.], subjects)
    counts = depths[:, None] * proportions
    data = ECGLMMData((counts,), (np.eye(paths),), np.ones((len(counts), 1)), clusters)
    return data, labels, clusters, proportions, depths


@pytest.mark.parametrize("covariance_source", ["likelihood", "posterior"])
def test_identity_mapping_recovers_effective_depth(covariance_source):
    data, labels, subjects, proportions, depths = direct_data(paths=3)
    result = quantify_effective_paths(data, [0, 1, 2], labels, subjects, baseline=np.ones(3) / 3, covariance_source=covariance_source)
    assert np.allclose(result["effective_depths"], depths, rtol=1e-7)
    expected = (depths[:, None] * proportions + .25 / 3) / (depths[:, None] + .25)
    assert np.allclose(result["proportions"], expected, atol=2e-5)
    assert np.allclose(result["counts"].sum(axis=1), depths)


def test_joint_dm_null_subject_offsets_are_not_a_cell_type_effect():
    data, labels, subjects, _, _ = direct_data()
    result = quantify_effective_paths(data, [0, 1], labels, subjects, baseline=np.ones(2) / 2)
    tested = joint_path_dm_test(result)
    assert tested["converged"]
    assert tested["p_value"] > .5
    assert tested["n_subjects"] == 12
    assert tested["standardized_means"].shape == (2, 2)


def test_joint_dm_detects_effect_and_reverses_standardized_direction():
    data, labels, subjects, _, _ = direct_data(subjects=20)
    counts = data.counts[0].copy()
    totals = counts.sum(axis=1)
    counts[labels == 0] = totals[labels == 0, None] * [.7, .3]
    counts[labels == 1] = totals[labels == 1, None] * [.3, .7]
    data = ECGLMMData((counts,), data.compatibility, data.design, data.clusters)
    quantified = quantify_effective_paths(data, [0, 1], labels, subjects, baseline=np.ones(2) / 2)
    result = joint_path_dm_test(quantified)
    assert result["converged"] and result["p_value"] < 1e-6
    reverse = joint_path_dm_test(quantified, labels=1 - labels)
    assert reverse["converged"] and reverse["p_value"] < 1e-6
    assert np.allclose(result["standardized_means"], reverse["standardized_means"][::-1], atol=1e-5)


def test_invalid_covariance_source_is_rejected():
    data, labels, subjects, _, _ = direct_data()
    with pytest.raises(ValueError, match="covariance source"):
        quantify_effective_paths(data, [0, 1], labels, subjects, covariance_source="pretend_counts")


def test_failed_joint_dm_variant_keeps_null_and_observed_denominators(monkeypatch):
    from extra_scripts import audit_path_reporting_omnibus as audit
    data, labels, subjects, _, _ = direct_data()
    def failed(*args, **kwargs):
        raise ValueError("test failure")
    monkeypatch.setattr(audit, "quantify_effective_paths", failed)
    observed, null, details = audit.joint_dm_reports(data, np.ones(2) / 2, np.array([0, 1]), labels, subjects, 8, 0)
    assert len(observed) == 3
    assert len(null) == 24
    assert all(result["p_value"] == 1 and not result["converged"] for result in observed.values())
    assert np.isnan(details["adjusted_effects"]).all()


def test_joint_paired_assessment_preserves_model_direction():
    import json
    import pandas as pd
    from extra_scripts.assess_joint_path_dm import paired_table
    frame = pd.DataFrame([{"test_id": "b_a", "block_id": "b", "gene_id": "gene.1", "level_a": "A", "level_b": "B", "levels": json.dumps(["A", "B"]), "path_signatures": json.dumps([[1], [2]]), "p_value": .01, "fdr": .03, "converged": True, "median_gene_umis": 50., "adjusted_effects": json.dumps([[0, 0], [-.4, .4]]), "standardized_means": json.dumps([[.7, .3], [.3, .7]])}])
    result = paired_table(frame)
    assert result.loc[0, "gene_id"] == "gene"
    assert result.loc[0, "pair_id"] == "A||B"
    assert result.loc[0, "effect_vector"] == [-.4, .4]


def test_joint_statistical_shards_do_not_invoke_reporting_fallback(tmp_path):
    import pandas as pd
    from extra_scripts.summarize_path_reporting_omnibus import load_shards
    root = tmp_path / "pairwise_fold0"
    (root / "shard_0").mkdir(parents=True)
    pd.DataFrame([{"strategy": "joint DM LRT A0.25 likelihood", "test_id": "b_a", "p_value": .05}]).to_csv(root / "shard_0/observed.tsv", sep="\t", index=False)
    result = load_shards(root, expected=1)
    assert len(result) == 1
    assert "effect" not in result


def test_multilevel_joint_dm_tests_all_cell_types_and_is_path_equivariant():
    from scipy.special import softmax
    labels = np.tile(np.arange(3), 12)
    subjects = np.repeat(np.arange(12), 3)
    rng = np.random.default_rng(37)
    subject_logits = rng.normal(0, .25, (12, 3))
    logits = subject_logits[subjects] + np.asarray([[.8, 0, -.8], [0, .8, -.8], [-.8, 0, .8]])[labels]
    counts = 300 * softmax(logits, axis=1)
    quantified = {"counts": counts, "subjects": subjects, "labels": labels, "effective_depths": np.full(36, 300.)}
    result = joint_path_dm_test(quantified)
    assert result["converged"] and result["p_value"] < 1e-8
    assert result["degrees_of_freedom"] == 4
    assert result["standardized_means"].shape == (3, 3)
    assert np.allclose(result["standardized_means"].sum(axis=1), 1)
    order = [2, 0, 1]
    reverse = joint_path_dm_test({**quantified, "counts": counts[:, order]})
    assert reverse["converged"]
    assert np.allclose(reverse["standardized_means"], result["standardized_means"][:, order], atol=2e-5)
    assert reverse["statistic"] == pytest.approx(result["statistic"], rel=1e-5)


def test_failed_omnibus_fit_cannot_retain_an_extreme_p_value():
    import json
    import pandas as pd
    from extra_scripts.assess_omnibus_inference_audit import complete_failed_tests
    reference = pd.DataFrame([{"block_id": "b", "test_id": "b", "adjusted_effects": json.dumps([[0, 0], [.4, -.4]])}])
    failed = reference.assign(strategy="joint DM", converged=False, p_value=1e-10, raw_p_value=1e-12, fdr=1e-9, statistic=50.)
    result = complete_failed_tests(failed, reference)
    assert not result.loc[0, "fit_available"]
    assert result.loc[0, "p_value"] == 1
    assert result.loc[0, "raw_p_value"] == 1
    assert result.loc[0, "fdr"] == 1
    assert result.loc[0, "statistic"] == 0
    assert np.isnan(np.asarray(json.loads(result.loc[0, "adjusted_effects"]))).all()
