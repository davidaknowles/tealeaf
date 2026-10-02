"""Check complete event mapping and the hybrid's inclusion-effect convention."""

import sys

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts import assess_event_tilgner_replication as audit
from extra_scripts import retest_suppa2_hybrid_psi
from extra_scripts.publish_suppa2_hybrid_lr_fix import reporting_table
from extra_scripts.plot_tilgner_method_replication import rank_agreement_table
from tealeaf.sc import differential


def test_hybrid_class_likelihood_and_total_pseudocount():
    baseline = np.array([0.45, 0.45, 0.10])
    counts = (np.array([80., 20., 10.]), np.array([40., 10., 5.]))
    fit = differential.fit_path_perturbation(
        counts, (np.eye(3), np.eye(3)), baseline, np.array([0, 1, -1]),
        path_pseudocount=32, path_pseudocount_scaling="total",
    )
    expected = (120 + 16) / (150 + 32)
    assert fit.converged
    np.testing.assert_allclose(fit.path_proportions, [expected, 1 - expected], atol=2e-5)
    np.testing.assert_allclose(fit.theta.sum(), 1)
    np.testing.assert_allclose(fit.theta[:2].sum(), 0.90)
    np.testing.assert_allclose(fit.theta[2], 0.10)
    np.testing.assert_allclose(fit.path_logratios[0], np.log(expected / (1 - expected)) / np.sqrt(2), atol=1e-4)


def test_all_tests_mapped_without_native_selection(tmp_path, monkeypatch):
    # The source has higher inclusion in B in both biological replicates.
    source = sparse.csr_matrix([[10., 10., 80., 80.], [90., 90., 20., 20.]])
    features = pd.DataFrame({"transcript_id": ["inc", "exc"]})
    columns = pd.DataFrame({"tealeaf_cell_type": ["A", "A", "B", "B"], "replicate": [1, 2, 1, 2], "column": range(4)})
    monkeypatch.setattr(audit, "read_tilgner_matrix", lambda *args: (source, features, columns))
    events = tmp_path / "events.tsv"
    pd.DataFrame({"feature_id": ["SUPPA2:event"], "included": ["inc"], "excluded": ["exc"]}).to_csv(events, sep="\t", index=False)
    tests = tmp_path / "tests.tsv"
    records = []
    for index in range(3):
        records.append({"method": "Tealeaf EC; SUPPA2 event definitions", "effect": "cell_type", "contrast_id": str(index), "level_a": "A", "level_b": "B", "feature_id": "SUPPA2:event", "event_type": "SE", "p_value": .01, "raw_p_value": .001 + index / 10000, "statistic": 20 - index, "effect_size": 1 if index != 1 else -1, "converged": index != 2})
    pd.DataFrame(records).to_csv(tests, sep="\t", index=False)
    output, summary = tmp_path / "audit.tsv.gz", tmp_path / "summary.tsv"
    monkeypatch.setattr(sys, "argv", ["audit", "--tests", str(tests), "--events", str(events), "--tilgner-matrix", str(tmp_path), "--gtf", str(tmp_path / "unused.gtf"), "--output", str(output), "--summary", str(summary), "--top-per-contrast", "0"])
    audit.main()
    mapped = pd.read_csv(output, sep="\t")
    assert len(mapped) == 2
    assert mapped.pooled_replicated.tolist() == [True, False]
    np.testing.assert_allclose(mapped.long_read_effect, .70)
    np.testing.assert_allclose(mapped.raw_p_value, [.001, .0011])
    np.testing.assert_allclose(mapped.statistic, [20, 19])
    assert pd.read_csv(summary, sep="\t").scope.eq("all converged tests").all()


def test_profiled_event_mass_prevents_nuisance_change_reversing_inclusion():
    mapping = np.array([[1., 0., 0.], [0., 1., .8], [0., 0., .2]])
    baseline = np.array([.4, .4, .2])
    true_a, true_b = np.array([.48, .32, .2]), np.array([.16, .04, .8])
    fixed, profiled = [], []
    for theta in (true_a, true_b):
        counts = (10000 * (mapping @ theta),)
        fixed.append(differential.fit_path_perturbation(counts, (mapping,), baseline, np.array([0, 1, -1])).path_proportions[0])
        fit = differential.fit_event_path_perturbation(counts, (mapping,), baseline, np.array([0, 1, -1]))
        assert fit.converged
        np.testing.assert_allclose(fit.theta, theta, atol=2e-4)
        profiled.append(fit.path_proportions[0])
    assert fixed[1] < fixed[0]
    np.testing.assert_allclose(profiled[1] - profiled[0], .2, atol=2e-4)


def test_profiled_mass_and_inclusion_prior_have_separate_targets():
    counts = (np.array([80., 20., 100.]), np.array([40., 10., 50.]))
    fit = differential.fit_event_path_perturbation(counts, (np.eye(3), np.eye(3)), np.array([.45, .45, .1]), np.array([0, 1, -1]), path_pseudocount=32, path_pseudocount_scaling="total")
    expected = (120 + 16) / (150 + 32)
    assert fit.converged
    np.testing.assert_allclose(fit.path_proportions[0], expected, atol=2e-5)
    np.testing.assert_allclose(fit.theta[:2].sum(), .5, atol=2e-5)
    assert fit.covariance.identifiable
    np.testing.assert_allclose(fit.covariance.covariance[0, 0], 1 / (2 * 182 * expected * (1 - expected)), rtol=2e-4)


def test_event_inclusion_fit_is_not_limited_by_extreme_baseline():
    fit = differential.fit_event_path_perturbation((np.array([800., 200.]),), (np.eye(2),), np.array([1e-12, 1 - 1e-12]), np.array([0, 1]), path_pseudocount=64, path_pseudocount_scaling="total")
    assert fit.converged
    np.testing.assert_allclose(fit.path_proportions[0], (800 + 32) / (1000 + 64), atol=2e-5)


def test_transcript_information_is_finite_for_tiny_ec_masses():
    theta = np.array([.4, .6])
    information = differential.transcript_fisher_information(theta, (np.eye(2) * 1e-310,), (100.,))
    assert np.isfinite(information).all()
    np.testing.assert_allclose(information, 100 * (np.diag(theta) - np.outer(theta, theta)), rtol=1e-6)


def test_psi_sensitivity_pairs_saved_numeric_celltype_labels(tmp_path, monkeypatch):
    shard = tmp_path / "shards" / "shard_0"
    shard.mkdir(parents=True)
    record = {"test_id": "event", "block_id": "block", "level_a": "A", "level_b": "B", "report_pseudocount": 1}
    pd.DataFrame([record]).to_csv(shard / "paired_path.tsv", sep="\t", index=False)
    usage = [{"test_id": "event", "subject": subject, "cell_type": level, "inclusion": .2 if level == 0 else .5 + subject / 100} for subject in range(4) for level in (0, 1)]
    pd.DataFrame(usage).to_csv(shard / "path_usage.tsv.gz", sep="\t", index=False)
    output = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", ["retest", "--shards", str(shard.parent), "--output-dir", str(output), "--null-replicates", "2"])
    retest_suppa2_hybrid_psi.main()
    result = pd.read_csv(output / "paired_path.tsv", sep="\t").iloc[0]
    assert result.converged and result.n_subjects == 4
    np.testing.assert_allclose(result.effect_size, .315)
    assert result.p_value < .001
    assert len(pd.read_csv(output / "paired_path_null.tsv.gz", sep="\t")) == 2


def test_reporting_effect_does_not_change_testing_fields():
    original = pd.DataFrame({"effect_size": [.2], "report_psi_effect": [-.05], "p_value": [.001], "mean_difference_norm": [.2], "statistic": [30.]})
    reported = reporting_table(original)
    assert reported.effect_size.iloc[0] == -.05
    assert reported.test_ilr_effect_size.iloc[0] == .2
    assert reporting_table(reported).test_ilr_effect_size.iloc[0] == .2
    pd.testing.assert_frame_equal(reported[["p_value", "statistic", "mean_difference_norm"]], original[["p_value", "statistic", "mean_difference_norm"]])


def test_hybrid_rank_uses_own_complete_universe_and_continuous_ties(tmp_path):
    events = pd.DataFrame({"method": ["Tealeaf EC; SUPPA2 event definitions"] * 250, "feature_id": [f"event_{i}" for i in range(250)], "p_value": .01, "raw_p_value": np.arange(250, 0, -1) / 1e6, "statistic": np.arange(250), "pooled_replicated": True, "mapping_complete": True, "minimum_pooled_depth": 100})
    base, mapped = tmp_path / "base.tsv", tmp_path / "mapped.tsv"
    events.iloc[:0].to_csv(base, sep="\t", index=False)
    events.to_csv(mapped, sep="\t", index=False)
    ranked = rank_agreement_table(base, None, None, event_replication_paths=[mapped])
    assert len(ranked) == 200
    assert ranked.method.eq("Tealeaf/SUPPA2 hybrid").all()
    assert ranked.feature_id.iloc[0] == "event_249"
    assert ranked.feature_id.iloc[-1] == "event_50"
    assert ranked["rank"].tolist() == list(range(1, 201))
