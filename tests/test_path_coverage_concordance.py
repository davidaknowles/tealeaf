import pandas as pd

from extra_scripts.plot_microglia_path_coverage_concordance import make_feature_contrasts, summarize_agreement, summarize_block_contrasts, summarize_rank_concordance


def test_feature_contrast_uses_common_subjects_and_support_thresholds():
    rows = []
    for subject, low, high in (("s1", 0.1, 0.7), ("s2", 0.2, 0.8), ("s3", 0.3, 0.9)):
        for cell_type, fitted, observed in (("low", 0.2, low), ("high", 0.8, high)):
            rows.append({"test_id": "b1", "subject": subject, "cell_type": cell_type, "primer": "poly(dT)", "feature_id": "E:1-2", "feature_type": "exon", "start": 1, "end": 2, "path_numbers": "1", "fitted_inclusion": fitted, "observed_inclusion": observed, "evidence_denominator": 4})
    contrasts = make_feature_contrasts(pd.DataFrame(rows))
    assert len(contrasts) == 1
    assert contrasts.loc[0, "n_common_subjects"] == 3
    assert abs(contrasts.loc[0, "fitted_delta"] - 0.6) < 1e-12
    assert abs(contrasts.loc[0, "observed_delta"] - 0.6) < 1e-12


def test_block_and_agreement_summaries():
    contrasts = pd.DataFrame([
        {"test_id": test_id, "feature_type": "exon", "primer": "poly(dT)", "feature_id": feature_id, "fitted_delta": fitted, "observed_delta": observed}
        for test_id, fitted, observed in (("b1", 0.2, 0.1), ("b2", 0.4, 0.3), ("b3", 0.8, 0.7))
        for feature_id in ("E1", "E2")
    ])
    blocks = summarize_block_contrasts(contrasts)
    assert len(blocks) == 3
    assert blocks["n_features"].eq(2).all()
    summary = summarize_agreement(blocks, bootstrap_replicates=20)
    assert summary.loc[0, "n_blocks"] == 3
    assert summary.loc[0, "sign_agreement"] == 1
    assert summary.loc[0, "spearman"] == 1


def test_rank_concordance_pools_raw_numerators_and_denominators():
    rows = []
    for cell_type, fitted, observed in (("A", 0.1, 0.2), ("B", 0.5, 0.6), ("C", 0.9, 0.8)):
        for subject in ("s1", "s2"):
            rows.append({"test_id": "b1", "feature_id": "E:1-2", "feature_type": "exon", "primer": "poly(dT)", "cell_type": cell_type, "subject": subject, "fitted_inclusion": fitted, "observed_inclusion": observed, "unclipped_observed_ratio": observed, "evidence_denominator": 10})
    features, blocks, summary = summarize_rank_concordance(pd.DataFrame(rows), bootstrap_replicates=20)
    assert features.loc[0, "rank_correlation"] == 1
    assert blocks.loc[0, "rank_correlation"] == 1
    assert summary.loc[0, "positive_fraction"] == 1
