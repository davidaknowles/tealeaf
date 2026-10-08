import json
import sys

import numpy as np
import pandas as pd
import pytest

from extra_scripts import audit_event_score_subject_influence
from tealeaf.sc.path_score_mixed import MODEL_VERSION, binary_score_components_from_records, mixed_score_test


def archive_fixture(tmp_path, missing=False):
    shard = tmp_path / "cohort" / "shard_0"
    shard.mkdir(parents=True)
    records, observed, null = [], [], []
    for test_id, information, mean in (("a", [1e6, 1., 1., 1., 1.], .5), ("b", [4.] * 5, .25)):
        local = [dict(test_id=test_id, subject=f"s{index}", score=value * mean, information=value, biological_shape=4., reference_information=value, report_inclusion_a=.2, report_inclusion_b=.3) for index, value in enumerate(information)]
        components = binary_score_components_from_records(local)
        fitted = mixed_score_test(components.scores, components.information, components.biological_shapes, reference_information=components.reference_information, scalar_fast=True)
        observed.append(dict(test_id=test_id, gene_id=f"g{test_id}", n_subjects=5, n_isoforms=2, n_ecs=3, median_gene_umis=100., p_value=fitted["p_value"]))
        if not missing or test_id != "a":
            records.extend(local)
        null.extend(dict(test_id=test_id, replicate=replicate, p_value=fitted["p_value"] if test_id == "a" else .5) for replicate in range(32))
    pd.DataFrame(records).to_csv(shard / "subject_scores.tsv.gz", sep="\t", index=False)
    pd.DataFrame(observed).to_csv(shard / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame(null).to_csv(shard / "paired_path_null.tsv.gz", sep="\t", index=False)
    (shard / "summary.json").write_text(json.dumps(dict(completed=2, failures=0, tests_in_shard=2)))
    (shard / "settings.json").write_text(json.dumps(dict(model_version=MODEL_VERSION, arguments=dict(information_metric="reference"))))
    merged = tmp_path / "merged"
    merged.mkdir()
    table = pd.DataFrame(observed).rename(columns={"p_value": "raw_p_value"})
    table["p_value"], table["complete_subject_fits"] = [.004, .5], True
    table.to_csv(merged / "paired_path.tsv", sep="\t", index=False)
    return shard.parent, merged


@pytest.mark.parametrize("missing", [False, True])
@pytest.mark.parametrize("compare", [False, True])
def test_complete_archive_diagnostic_retains_every_request(tmp_path, monkeypatch, missing, compare):
    cohort, merged = archive_fixture(tmp_path, missing)
    output = tmp_path / "output"
    argv = ["audit", "--cohort-root", str(cohort), "--merged-dir", str(merged), "--output-dir", str(output), "--shard-count", "1"]
    monkeypatch.setattr(sys, "argv", argv + (["--compare-proportions"] if compare else []))
    audit_event_score_subject_influence.main()
    cases = pd.read_csv(output / "diagnostics.tsv.gz", sep="\t").set_index("test_id")
    assert set(cases.index) == {"a", "b"}
    assert cases.at["b", "effective_weighted_subjects"] == pytest.approx(5.)
    if compare:
        assert cases.at["b", "proportion_status"] == "ok"
        assert cases.at["b", "proportion_fitted_mean"] == pytest.approx(.125)
        assert cases.at["b", "proportion_fitted_p_value"] == pytest.approx(cases.at["b", "fitted_p_value"], rel=1e-7)
    if missing:
        assert cases.at["a", "status"] == "missing_archive"
    else:
        assert cases.at["a", "status"] == "ok"
        assert cases.at["a", "maximum_precision_share"] > .9999
        assert cases.at["a", "n_sign_draws_below_threshold"] == 32
        assert cases.at["a", "leave_dominant_out_available"]
        assert cases.at["a", "leave_dominant_out_p_value"] > .05
    receipt = json.loads((output / "manifest.json").read_text())
    assert receipt["requested"] == 2 and receipt["declared_cohort_tests"] == 2
    assert "not LR" in receipt["selection"]
    if compare:
        summary = pd.read_csv(output / "summary.tsv", sep="\t").set_index("panel")
        assert summary.at["small observed raw p-value", "proportion_valid"] == (0 if missing else 1)
        assert summary.at["subject-count-matched weak control", "proportion_valid"] == 1
        assert summary.proportion_mean_outside_simplex_difference.sum() == 0
