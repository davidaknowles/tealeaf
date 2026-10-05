import pandas as pd
import pytest
from pathlib import Path

from tealeaf.sc.junction_benchmark import index_subject_paired_contrasts


def example():
    samples = pd.DataFrame({"sample_id": ["a1", "a2", "b1", "b2"], "subject": ["s1", "s2", "s1", "s2"], "cell_type": ["a", "a", "b", "b"]})
    contrast = {"level_a": "a", "level_b": "b", "samples_a": ["a2", "a1"], "samples_b": ["b1", "b2"], "paired_subjects": ["s2", "s1"]}
    mapping = {"a1": 10, "a2": 20, "b1": 30, "b2": 40}
    return samples, contrast, mapping


def test_count_indices_follow_subject_pairing_not_manifest_row_order():
    samples, contrast, mapping = example()
    result, = index_subject_paired_contrasts([contrast], samples, mapping)
    assert result["paired_subjects"] == ["s1", "s2"]
    assert result["indices_a"] == [10, 20]
    assert result["indices_b"] == [30, 40]
    assert contrast["samples_a"] == ["a2", "a1"]


def test_unpaired_or_mislabeled_samples_are_rejected():
    samples, contrast, mapping = example()
    samples.loc[3, "subject"] = "s3"
    with pytest.raises(ValueError, match="different subjects"):
        index_subject_paired_contrasts([contrast], samples, mapping)
    samples.loc[3, "subject"] = "s2"
    samples.loc[3, "cell_type"] = "a"
    with pytest.raises(ValueError, match="cell types"):
        index_subject_paired_contrasts([contrast], samples, mapping)


def test_duplicate_indices_or_subjects_are_rejected():
    samples, contrast, mapping = example()
    mapping["b2"] = 30
    with pytest.raises(ValueError, match="unique nonnegative"):
        index_subject_paired_contrasts([contrast], samples, mapping)
    mapping["b2"] = 40
    samples.loc[1, "subject"] = "s1"
    with pytest.raises(ValueError, match="multiple pseudobulks"):
        index_subject_paired_contrasts([contrast], samples, mapping)


def test_rmats_direction_loading_converts_native_delta_and_keeps_all_tests(tmp_path):
    import numpy as np
    from extra_scripts.audit_split_coverage_direction import load_paired_rmats

    table = pd.DataFrame({"method": ["rMATS"] * 4, "contrast_id": ["a_b"] * 4, "effect": ["cell_type"] * 4, "stratum": ["all"] * 4, "level_a": ["a"] * 4, "level_b": ["b"] * 4, "feature_id": ["e1", "e2", "e3", "e4"], "p_value": [.001, .8, .0001, np.nan], "q_value": [.05, 1, .01, np.nan], "effect_size": [.3, -.2, .8, .5], "gene_id": ["g.2", "g.2", "h", "g"], "statistical_model": ["PAIRADISE"] * 4, "effect_orientation": ["a_minus_b"] * 4})
    path = tmp_path / "tests.tsv.gz"
    table.to_csv(path, sep="\t", index=False)
    result = load_paired_rmats(path, pd.DataFrame({"gene_id": ["g"]}))
    assert result.feature_id.tolist() == ["e1", "e2"]
    assert np.allclose(result.effect_size, [-.3, .2])
    assert result.effect_orientation.eq("b_minus_a").all()
    table["statistical_model"] = "unpaired"
    table.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="native paired"):
        load_paired_rmats(path)


def test_rmats_retries_preserve_existing_work(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import extra_scripts.run_rmats_comparison as runner

    previous = tmp_path / "fold0_0"
    previous.mkdir()
    marker = previous / "marker.txt"
    marker.write_text("previous work")
    monkeypatch.setattr(runner, "run_command", lambda *args, **kwargs: None)
    monkeypatch.setattr(runner, "parse_result", lambda *args: pd.DataFrame({"p_value": [.5]}))
    contrast = {"index": 0, "indices_a": [0], "indices_b": [1], "paired_subjects": ["s"]}
    args = SimpleNamespace(paired_stats=True, paired_jcec_only=False)
    result = runner.process_one((0, contrast, tmp_path, tmp_path, 1, args))
    assert result[-1] is None
    assert marker.read_text() == "previous work"
    saved = pd.read_csv(result[3], sep="\t")
    assert saved.statistical_model.iloc[0] == "PAIRADISE"
    assert saved.effect_orientation.iloc[0] == "a_minus_b"


def test_paired_output_join_rejects_row_count_mismatch(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import extra_scripts.run_rmats_comparison as runner

    (tmp_path / "JCEC.raw.input.SE.txt").write_text("counts")
    (tmp_path / "fromGTF.SE.txt").write_text("events")

    def command_stub(command, *args, **kwargs):
        if str(runner.INCLUSION_LEVEL) in map(str, command):
            Path(command[-1]).write_text("inclusion_header\n")
        if str(runner.PAIRED_MODEL) in map(str, command):
            Path(command[-1]).write_text("pvalue_header\nextra_row\n")

    monkeypatch.setattr(runner, "run_command", command_stub)
    with pytest.raises(ValueError, match="different row counts"):
        runner.run_paired_jcec_only(tmp_path, SimpleNamespace(rscript=None), 1)


def test_direction_table_does_not_report_empty_selection_as_zero_agreement():
    from extra_scripts.publish_rmats_split_data import direction_cell, render_table1_row

    assert direction_cell(pd.Series({"n_events": 0})) == "--- (0 events)"
    assert direction_cell(pd.Series({"n_events": 5, "n_direction_evaluable": 0})) == "--- (unevaluable)"
    metrics = pd.DataFrame({"method": ["Tealeaf", "rMATS (paired)"], "shared_gene_pairs": [10, 11], "shared_genes": [4, 4]}).set_index("method")
    with pytest.raises(ValueError, match="shared universes differ"):
        render_table1_row(metrics, None)


def test_publisher_updates_all_historical_claims_and_guards_user_edits(tmp_path, monkeypatch):
    import hashlib
    import extra_scripts.publish_rmats_split_data as publisher

    audit = tmp_path / "analyses/split_coverage_direction"
    output = tmp_path / "analyses/comparator_suppa_rmats/rmats"
    docs = tmp_path / "docs"
    for directory in (audit, output, docs):
        directory.mkdir(parents=True)
    (audit / "README.md").write_text("The historical rMATS row is retained separately, its paired full-data rMATS results are not a substitute for split effects.\n")
    notebook = tmp_path / "LABNOTEBOOK.md"
    notebook.write_text("Final split metrics are pending, the archived full-data long-read analysis remains unchanged.\n")
    original = "\n".join([
        r"rMATS (native, historical) & old statistics\\",
        "Asterisks mark unavailable historical rMATS split event tables, not zero agreement; full-data effects are not substituted.",
        "The rMATS row is explicitly the historical unpaired split audit. The corrected paired-JCEC rMATS rerun is integrated below as a full-data event comparator, rather than being presented as a paired split result that was not rerun.",
        "The historical rMATS split event tables are unavailable, so its aggregate row is retained without fabricated directional values or substituted full-data effects.",
        "The split reproducibility row for rMATS above is retained as a historical unpaired audit; a new split paired rerun would be required to replace it with a directly matched conjunction analysis.",
        "The new paired-JCEC rMATS run restores 1,497 full-data event FDR calls and is included in the event-level table and long-read comparison above, but those calls are not substituted for the historical split row because the paired split manifests were not rerun.",
        r"scQuint & \(-.425 / -.338\) & \(-.069 / -.080\)\\",
    ]).encode()
    document = docs / "differential.tex"
    document.write_bytes(original)
    digest = hashlib.sha256(original).hexdigest()
    methods = ["Tealeaf", "rMATS (paired)"]
    metrics = pd.DataFrame({"comparison": ["rMATS (paired)"] * 2, "method": methods, "shared_gene_pairs": [2000] * 2, "shared_genes": [400] * 2, "replicated_bh": [100, 20], "heldout_nominal_replication": [.9, .8], "spearman_logp": [.8, .6]})
    metrics.to_csv(audit / "table1_gene_metrics.tsv", sep="\t", index=False)
    direction = pd.DataFrame({"comparison": ["rMATS (paired)"] * 2, "method": methods, "selection": ["event_BH_union"] * 2, "coverage_bin": ["all"] * 2, "n_events": [100, 100], "n_direction_evaluable": [90, 100], "agree": [89, 95], "agreement": [89 / 90, .95]})
    direction.to_csv(audit / "direction_summary.tsv", sep="\t", index=False)
    pd.DataFrame({"comparison": ["rMATS (paired)"] * 4, "unit": ["gene"] * 4, "method": methods * 2, "fold": [0, 0, 1, 1], "rho_p_coverage": [-.2, -.3, -.4, -.5]}).to_csv(audit / "coverage_correlations.tsv", sep="\t", index=False)
    monkeypatch.setattr("sys.argv", ["publish", "--repo-root", str(tmp_path), "--expected-document-sha", "wrong"])
    with pytest.raises(ValueError, match="manuscript changed"):
        publisher.main()
    assert document.read_bytes() == original
    monkeypatch.setattr("sys.argv", ["publish", "--repo-root", str(tmp_path), "--expected-document-sha", digest])
    publisher.main()
    changed = document.read_text()
    assert "rMATS (paired JCEC) & 2,000 & 400" in changed
    assert "rMATS (native, historical)" not in changed
    assert "rMATS paired JCEC & \\(-0.200 / -0.400\\)" in changed
    assert "because the paired split manifests were not rerun" not in changed
    assert "split event tables are unavailable" not in changed
    assert "Completed all 380 paired split contrasts" in notebook.read_text()
    assert "now replace the historical unpaired row" in (audit / "README.md").read_text()
