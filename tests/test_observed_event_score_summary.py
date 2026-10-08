import json

import pandas as pd
import pytest

from extra_scripts.summarize_observed_event_score import collate_source, summarize, validate_source_identities


def test_observed_score_summary_preserves_failures_and_external_zeros(tmp_path):
    frame = pd.DataFrame(dict(scope=["native diagnostic"] * 3, test_id=["a", "b", "c"], source=["original_binary"] * 3, p_value=[.01, 1., .02], converged=[True, False, True], report_complete=[True, False, True], report_effect=[.2, None, .2], score_effect=[.2, None, .2], chi_square_p_value=[.001, None, .001], pooled_inclusion=[.5, .999, .5], pooled_event_mass=[.1, .001, .1], long_read_effect=[.3, .3, 0.]))
    settings = dict(shard_count=1, n_requested=3)
    shard = tmp_path / "shard_0"
    shard.mkdir()
    (shard / "settings.json").write_text(json.dumps(settings))
    frame.to_csv(shard / "tests.tsv", sep="\t", index=False)
    table, _ = collate_source(tmp_path)
    row = summarize(table).iloc[0]
    assert row.n_requested == 3 and row.n_complete_inference == 2
    assert row.n_report_agrees == row.n_score_agrees == 1
    assert row.n_external_nonzero == 2
    frame.loc[1, "p_value"] = .01
    frame.to_csv(shard / "tests.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="p1"):
        collate_source(tmp_path)


def test_source_specific_random_panels_require_explicit_flag_and_same_native_panel():
    def frame(prefix):
        return pd.DataFrame([dict(scope="fixed random real-data diagnostic", test_id=f"{prefix}_{index}") for index in range(32)] + [dict(scope="native top100 diagnostic, not an inference family", test_id=f"native_{index}") for index in range(100)])

    parts = [frame("first"), frame("second")]
    with pytest.raises(ValueError, match="same declared"):
        validate_source_identities(parts)
    assert not validate_source_identities(parts, True)
    assert validate_source_identities([parts[0], parts[0]])
    parts[1].loc[131, "test_id"] = "different native hit"
    with pytest.raises(ValueError, match="identical native"):
        validate_source_identities(parts, True)
