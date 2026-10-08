import hashlib
import json
import sys

import numpy as np
import pandas as pd
import pytest

from extra_scripts import audit_pooled_null_signal


def run_audit(tmp_path, monkeypatch, observed, null):
    observed_path, null_path = tmp_path / "observed.tsv", tmp_path / "null.tsv.gz"
    observed.to_csv(observed_path, sep="\t", index=False)
    null.to_csv(null_path, sep="\t", index=False)
    hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (observed_path, null_path)}
    output = tmp_path / "audit"
    monkeypatch.setattr(sys, "argv", ["audit", "--observed", str(observed_path), "--null", str(null_path), "--output-dir", str(output)])
    audit_pooled_null_signal.main()
    assert hashes == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in hashes}
    return output


def example_tables():
    observed = pd.DataFrame(dict(test_id=["a", "b"], calibration_stratum=["same", "same"], raw_p_value=[1e-6, .3], p_value=[.004, .45]))
    null = pd.DataFrame(dict(test_id=["a", "a", "b", "b"], calibration_stratum=["same"] * 4, raw_p_value=[1e-7, .9, .05, .8], replicate=[0, 1, 0, 1]))
    return observed, null


def test_tail_audit_preserves_inputs_and_labels_parent_thresholds(tmp_path, monkeypatch):
    observed, null = example_tables()
    output = run_audit(tmp_path, monkeypatch, observed, null)
    result = pd.read_csv(output / "tail_contributions.tsv", sep="\t")
    row = result.loc[result.threshold.eq(.001)].iloc[0]
    assert row.observed_tests == 2 and row.training_draws == 4 and row.training_parents == 2
    assert row.raw_observed_calls == 1 and row.calibrated_observed_calls == 0
    assert row.raw_null_tail_draws == 1 and row.tail_fraction_from_small_observed_raw_parents == 1.
    assert result.loc[result.threshold.eq(.05), "tail_fraction_from_small_observed_raw_parents"].iloc[0] == .5
    receipt = json.loads((output / "manifest.json").read_text())
    assert "not proof" in receipt["interpretation"]
    assert "upstream" in receipt["scope"]


@pytest.mark.parametrize("defect", ["duplicate_observed", "duplicate_null", "foreign", "stratum", "negative", "nonfinite"])
def test_tail_audit_rejects_unaligned_or_invalid_inputs(tmp_path, monkeypatch, defect):
    observed, null = example_tables()
    if defect == "duplicate_observed":
        observed.loc[1, "test_id"] = "a"
    elif defect == "duplicate_null":
        null.loc[1, "replicate"] = 0
    elif defect == "foreign":
        null.loc[0, "test_id"] = "other"
    elif defect == "stratum":
        null.loc[0, "calibration_stratum"] = "wrong"
    elif defect == "negative":
        null.loc[0, "raw_p_value"] = -.1
    else:
        observed.loc[0, "p_value"] = np.nan
    with pytest.raises(ValueError):
        run_audit(tmp_path, monkeypatch, observed, null)
