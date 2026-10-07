import json

import numpy as np
import pandas as pd
import pytest

from extra_scripts.assess_event_input_controls import guard_completed_shards
from extra_scripts.audit_hybrid_input_count_null import trial_header
from tealeaf.sc.replication_audit import complete_cluster_fit


def test_informative_clusters_are_not_mistaken_for_successful_subject_fits():
    fitted = {"converged": True, "n_subjects": 4, "n_fitted_subjects": 6}
    assert complete_cluster_fit(fitted, 6)
    assert not complete_cluster_fit(fitted, 7)
    assert not complete_cluster_fit({**fitted, "n_subjects": 3}, 6)
    assert not complete_cluster_fit({"converged": True, "n_subjects": 4}, 6)


def test_count_null_header_has_event_identity_for_shared_summarizer():
    record = {"test_id": "SUPPA2:g;SE:event|cell_type|A|B", "feature_id": "SUPPA2:g;SE:event", "gene_id": "g"}
    assert trial_header(record) == {"test_id": record["test_id"], "block_id": "g;SE:event", "gene_id": "g", "n_paths": 2}


def write_shard(root, incomplete=False, failures=False):
    shard = root / "shard_0"
    shard.mkdir(parents=True)
    table = pd.DataFrame({"test_id": ["t1", "t2"], "n_samples": [10, 10], "n_subjects": [5, 4 if incomplete else 5], "report_n_subjects": [4, 5], "converged": [True, True], "path_pseudocount": [32, 32], "report_pseudocount": [1, 1], "profile_event_mass": [True, True], "p_value": [.01, .02], "statistic": [5, 4], "effect_size": [.2, .3], "report_psi_effect": [.2, .3], "report_ilr_effect": [.2, .3]})
    table.to_csv(shard / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame([{"test_id": key, "replicate": replicate, "p_value": .5} for key in table.test_id for replicate in range(32)]).to_csv(shard / "paired_path_null.tsv.gz", sep="\t", index=False)
    errors = [{"test_id": "SUPPA2:g;SE:event|cell_type|A|B", "error": "test"}] if failures else []
    (shard / "failures.json").write_text(json.dumps(errors))
    (shard / "summary.json").write_text(json.dumps({"completed": 2, "failures": len(errors), "tests_in_shard": 2 + len(errors)}))
    return shard


def test_guards_preserve_failed_hypotheses_and_separate_reporting_failure(tmp_path):
    write_shard(tmp_path / "raw", incomplete=True, failures=True)
    checks = guard_completed_shards(tmp_path / "raw", tmp_path / "guarded", 32, shard_count=1)
    table = pd.read_csv(tmp_path / "guarded/shard_0/paired_path.tsv", sep="\t").set_index("test_id")
    null = pd.read_csv(tmp_path / "guarded/shard_0/paired_path_null.tsv.gz", sep="\t")
    assert len(table) == 3
    assert np.isnan(table.loc["t1", "effect_size"])
    assert table.loc["t1", "p_value"] == .01
    assert table.loc["t2", "p_value"] == 1
    assert table.loc["SUPPA2:g;SE:event|cell_type|A|B", "p_value"] == 1
    assert set(null.test_id) == {"t1"}
    assert checks[0]["complete_subject_fits"] == 1


def test_guard_rejects_incomplete_arrays_or_null_families(tmp_path):
    shard = write_shard(tmp_path / "raw")
    with pytest.raises(FileNotFoundError):
        guard_completed_shards(tmp_path / "raw", tmp_path / "out", 32, shard_count=2)
    null = pd.read_csv(shard / "paired_path_null.tsv.gz", sep="\t").iloc[:-1]
    null.to_csv(shard / "paired_path_null.tsv.gz", sep="\t", index=False)
    with pytest.raises(ValueError, match="null realizations"):
        guard_completed_shards(tmp_path / "raw", tmp_path / "out2", 32, shard_count=1)
