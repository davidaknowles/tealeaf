import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from extra_scripts.audit_shared_score_variance import prepare_panel
from extra_scripts.audit_shared_variance_count_null import assess, count_signed_values
from tealeaf.sc.path_score_mixed import MODEL_VERSION


def fixture(tmp_path):
    source = tmp_path / "original" / "shard_0"
    source.mkdir(parents=True)
    observations, subjects = [], []
    ids = [f"g{index}" for index in range(5)] + ["failed"]
    for draw in range(2):
        for key in ids:
            valid = key != "failed"
            observations.append(dict(test_id=key, gene_id=key, draw=draw, p_value=.5 if valid else 1., complete_subject_fits=valid, n_subjects=5 if valid else 0, n_expected_subjects=6 if valid else 0, n_fitted_subjects=6 if valid else 0))
            if valid:
                for index in range(6):
                    information = 4. if index != 2 else 1e-14
                    subjects.append(dict(test_id=key, draw=draw, subject=f"s{index}", score=information * (.2 if index % 2 == 0 else -.1), information=information, reference_information=8., biological_shape=4.))
    observed = pd.DataFrame(observations)
    observed.to_csv(source / "observed.tsv", sep="\t", index=False)
    pd.DataFrame(subjects).to_csv(source / "subject_null_diagnostics.tsv.gz", sep="\t", index=False)
    assessed = tmp_path / "validated"
    assessed.mkdir()
    observed.to_csv(assessed / "trials.tsv.gz", sep="\t", index=False)
    recipe = dict(requested_ids=ids, draws=2, expected_strategies=["fixture"], score_coordinate="ilr", information_metric="reference", mixed_score_version=MODEL_VERSION)
    receipt = dict(source=str(source.parent), source_recipe=recipe, shards=[dict(shard=0, diagnostics_sha256=hashlib.sha256((source / "subject_null_diagnostics.tsv.gz").read_bytes()).hexdigest(), observed_sha256=hashlib.sha256((source / "observed.tsv").read_bytes()).hexdigest())])
    (assessed / "manifest.json").write_text(json.dumps(receipt))
    return assessed, source


def test_count_pilot_preserves_failures_refits_all_nulls_and_checks_source(tmp_path):
    assessed, source = fixture(tmp_path)
    args = SimpleNamespace(source_assessment=assessed, output_dir=tmp_path / "output")
    assess(args)
    result = pd.read_csv(args.output_dir / "trials.tsv.gz", sep="\t")
    assert len(result) == 12
    assert result.loc[result.test_id.eq("failed"), ["p_value", "raw_p_value"]].eq(1.).all().all()
    variances = pd.read_csv(args.output_dir / "variances.tsv", sep="\t")
    assert len(variances) == 66
    assert variances.available_training_genes.eq(5).all()
    null = pd.read_csv(args.output_dir / "null.tsv.gz", sep="\t")
    assert len(null) == 5 * 2 * 32 and not null.duplicated(["test_id", "draw", "replicate"]).any()
    assert null.groupby(["test_id", "draw"]).replicate.apply(lambda value: set(value) == set(range(32))).all()
    (source / "observed.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="changed after source assessment"):
        assess(SimpleNamespace(source_assessment=assessed, output_dir=tmp_path / "different"))


def test_count_signs_reproduce_sequential_original_stream_before_rank_selection(tmp_path):
    import zlib
    _, source = fixture(tmp_path)
    records = pd.read_csv(source / "subject_null_diagnostics.tsv.gz", sep="\t")
    ids = ["g0", "g1"]
    groups = dict(tuple(records.loc[records.draw.eq(1) & records.test_id.isin(ids)].groupby("test_id", sort=False)))
    panel, table = prepare_panel(groups, ids)
    streams = [np.random.default_rng(381924 + zlib.crc32(key.encode()) + 1721) for key in ids]
    for replicate in range(32):
        expected = np.concatenate([stream.choice((-1., 1.), size=6) for stream in streams])
        np.testing.assert_array_equal(count_signed_values(panel, table, ids, 1, replicate), panel.values * expected[panel.source_positions])
