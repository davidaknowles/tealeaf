import numpy as np
import pandas as pd

from extra_scripts.assess_tilgner_long_read_replication import load_path_usage, block_feature_rows


def test_complete_reporting_guard_does_not_average_away_a_failed_subject(tmp_path):
    shard = tmp_path / "shard_0"
    shard.mkdir()
    table = pd.DataFrame({"test_id": ["a", "a", "b", "b"], "cell_type": ["C"] * 4, "path_number": [1] * 4, "proportion": [.2, np.nan, .3, .5]})
    table.to_csv(shard / "path_usage.tsv", sep="\t", index=False)
    permissive = load_path_usage(tmp_path, {"a", "b"}).set_index("test_id")
    strict = load_path_usage(tmp_path, {"a", "b"}, require_complete=True).set_index("test_id")
    assert permissive.loc["a", "proportion"] == .2
    assert np.isnan(strict.loc["a", "proportion"])
    assert strict.loc["b", "proportion"] == .4


def test_outside_transcript_is_not_wrapped_to_the_last_long_read_path():
    signatures = [[], [[3, 4]]]
    block = {"gene_id": "g.1", "transcripts": ["t0.1", "t1.2", "outside.1"], "path_index": [0, 1, -1], "path_signatures": signatures}
    features = pd.DataFrame({"stable_gene_id": ["g"] * 3, "transcript_id": ["t0", "t1", "outside"], "row": [0, 1, 2]})
    result = block_feature_rows(block, signatures, features)
    assert result.row.tolist() == [0, 1]
    assert result.path_number.tolist() == [1, 2]
