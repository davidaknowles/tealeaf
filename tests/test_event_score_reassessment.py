import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from extra_scripts.reassess_event_score_archive import reassess
from extra_scripts.run_suppa2_tealeaf_hybrid import mixed_event_record
from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_score_mixed import shared_path_score_components, binary_subject_score_records, binary_score_components_from_records, aggregate_path_scores, paired_score_reporting, MODEL_VERSION


def fixture():
    subjects = np.repeat(np.arange(6), 2)
    labels = np.tile([0, 1], 6)
    counts = np.tile([[30., 20., 40., 10.], [50., 20., 20., 10.]], (6, 1))
    base = ECGLMMData((counts,), (np.eye(4),), np.ones((12, 1)), subjects)
    paths, baseline = np.array([0, 0, 1, -1]), np.array([.4, .2, .3, .1])
    components = shared_path_score_components(base, paths, labels, subjects, baseline=baseline, reporting_concentration=1.)
    args = SimpleNamespace(inference="mixed-score", max_iter=300, null_multistart=False, report_pseudocount=1., null_replicates=2, seed=20260927, export_path_usage=True, information_metric="absolute")
    event = SimpleNamespace(event_id="gene;SE:example", feature_id="SUPPA2:gene;SE:example", event_type="SE")
    return base, paths, labels, subjects, baseline, components, args, event


def test_binary_archive_roundtrip_preserves_model_and_missing_reports():
    *_, components, args, event = fixture()
    records = binary_subject_score_records(components, "test")
    restored = binary_score_components_from_records(records)
    for metric in ("absolute", "reference"):
        original = aggregate_path_scores(components, information_metric=metric)
        replay = aggregate_path_scores(restored, information_metric=metric)
        np.testing.assert_allclose(replay["p_value"], original["p_value"], rtol=1e-12)
        np.testing.assert_allclose(replay["mean_difference"], original["mean_difference"], rtol=1e-12)
    records[0]["report_inclusion_a"] = np.nan
    missing = binary_score_components_from_records(records)
    assert not paired_score_reporting(missing)["complete"]
    assert len(missing.subject_ids) == len(components.subject_ids)
    with pytest.raises(ValueError, match="unique"):
        binary_score_components_from_records(records + records[:1])


def write_source(source):
    base, paths, labels, subjects, baseline, components, args, event = fixture()
    archive = ([], [])
    row, null, usage = mixed_event_record(base, paths, labels, subjects, baseline, event, "gene", ("A", "B"), np.full(12, 100), 4, args, components=components, score_archive=archive)
    source.mkdir()
    pd.DataFrame([row]).to_csv(source / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame(archive[0]).to_csv(source / "score_contexts.tsv.gz", sep="\t", index=False)
    pd.DataFrame(archive[1]).to_csv(source / "subject_scores.tsv.gz", sep="\t", index=False)
    (source / "settings.json").write_text(json.dumps(dict(arguments=vars(args), model_version=MODEL_VERSION)))
    failure = dict(test_id="SUPPA2:gene;SE:failed|cell_type|A|B", error="null fit failed")
    (source / "failures.json").write_text(json.dumps([failure]))
    (source / "summary.json").write_text(json.dumps(dict(tests_in_shard=2, completed=1, failures=1)))
    return components, args, row, failure


def test_complete_shard_replay_keeps_unfitted_failures_and_matches_reference(tmp_path):
    source, output = tmp_path / "source", tmp_path / "reference"
    components, args, original, failure = write_source(source)
    assert reassess(source, output, information_metric="reference") == "reused score archive"
    table = pd.read_csv(output / "paired_path.tsv", sep="\t")
    expected = aggregate_path_scores(components, information_metric="reference")
    assert table.iloc[0].test_id == original["test_id"]
    np.testing.assert_allclose(table.iloc[0].p_value, expected["p_value"], rtol=1e-10)
    np.testing.assert_allclose(table.iloc[0].effect_size, original["effect_size"], rtol=1e-12)
    assert json.loads((output / "failures.json").read_text()) == [failure]
    assert json.loads((output / "summary.json").read_text())["tests_in_shard"] == 2
    assert len(pd.read_csv(output / "paired_path_null.tsv.gz", sep="\t")) == args.null_replicates
    for name in ("score_contexts.tsv.gz", "subject_scores.tsv.gz"):
        assert (output / name).read_bytes() == (source / name).read_bytes()
    with pytest.raises(ValueError, match="new output"):
        reassess(source, output, information_metric="reference")


def test_replay_rejects_incomplete_subject_archive(tmp_path):
    source = tmp_path / "source"
    write_source(source)
    scores = pd.read_csv(source / "subject_scores.tsv.gz", sep="\t").iloc[:-1]
    scores.to_csv(source / "subject_scores.tsv.gz", sep="\t", index=False)
    with pytest.raises(ValueError, match="subject archive is incomplete"):
        reassess(source, tmp_path / "out", information_metric="reference")
