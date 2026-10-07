import json
import sys

import numpy as np
import pandas as pd
import pytest

from extra_scripts.summarize_binary_marginal_ec import main


def test_repeated_count_null_collation_keeps_failed_trials(tmp_path, monkeypatch):
    cache, output = tmp_path / "cache", tmp_path / "results"
    shard = cache / "shard_0"
    shard.mkdir(parents=True)
    scenarios = [{"scenario": name, "draw": draw} for draw in (0, 1) for name in ("common composition", "biological precision 20")]
    settings = {"selected_ids": ["a", "b"], "shard_count": 1, "requested_scenarios": scenarios}
    (shard / "settings.json").write_text(json.dumps(settings))
    rows = [{"test_id": test_id, **scenario, "converged": test_id == "a", "p_value": .001, "standardized_means": "[[.2,.8],[.8,.2]]"} for test_id in ("a", "b") for scenario in scenarios]
    table = pd.DataFrame(rows)
    table.to_csv(shard / "tests.tsv.gz", sep="\t", index=False)
    monkeypatch.setattr(sys, "argv", ["summary", "--cache", str(cache), "--output-dir", str(output)])
    main()
    result = pd.read_csv(output / "tests.tsv.gz", sep="\t")
    assert len(result) == 8
    assert result.loc[result.test_id.eq("b"), "p_value"].eq(1).all()
    assert result.loc[result.test_id.eq("b"), "standardized_means"].isna().all()
    summary = pd.read_csv(output / "summary.tsv", sep="\t")
    assert summary.n_requested.eq(4).all()
    assert summary.n_converged.eq(2).all()
    assert summary.native_reject_rate_0_05.eq(.5).all()
    assert np.isnan(summary.median_runtime_seconds).all()
    table.iloc[:-1].to_csv(shard / "tests.tsv.gz", sep="\t", index=False)
    with pytest.raises(ValueError, match="missing"):
        main()
