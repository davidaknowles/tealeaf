import json

import pandas as pd
import pytest

from extra_scripts.summarize_independent_block_nuisance import STRATEGIES, collate


def test_complete_draw_family_and_failure_policy(tmp_path):
    for shard in range(2):
        directory = tmp_path / f"shard_{shard}"
        directory.mkdir()
        settings = {"shard_count": 2, "draws": 4, "shard_index": shard, "output_dir": str(directory)}
        (directory / "settings.json").write_text(json.dumps(settings))
        rows = [{"draw": draw, "strategy": strategy, "true_delta": 0., "converged": draw != 0, "p_value": .001, "estimated_delta": .1, "runtime_seconds": 1.} for draw in range(shard, 4, 2) for strategy in STRATEGIES]
        pd.DataFrame(rows).to_csv(directory / "tests.tsv.gz", sep="\t", index=False)
    table, summary, _ = collate(tmp_path)
    assert len(table) == 12
    assert table.loc[table.draw.eq(0), "p_value"].eq(1).all()
    assert table.loc[table.draw.eq(0), "estimated_delta"].isna().all()
    assert summary.n_requested.eq(4).all()
    assert summary.n_converged.eq(3).all()
    assert summary.native_reject_rate_0_05.eq(.75).all()
    assert summary.direction_agreement.isna().all()
    path = tmp_path / "shard_1" / "tests.tsv.gz"
    pd.read_csv(path, sep="\t").iloc[:-1].to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing or duplicate"):
        collate(tmp_path)


def test_null_local_marginal_is_independent_of_other_block():
    import numpy as np

    # The A-origin read likelihood depends on A but not within-A B shares.
    mapping = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [3, 3, 0, 0]], dtype=float)
    a = .3
    probabilities = []
    for b in (.1, .9):
        weights = np.array([a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)])
        mass = mapping @ weights
        probabilities.append(mass / mass.sum())
    np.testing.assert_allclose(*probabilities, atol=1e-14)
