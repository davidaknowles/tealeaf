import numpy as np
import pandas as pd
import pytest

from extra_scripts.merge_paired_path_test import empirical_null_calibration, moderate_scalar_tests
from types import SimpleNamespace


def test_moderate_scalar_tests_refits_each_null_family():
    test_ids = [f"test_{index}" for index in range(6)]
    table = pd.DataFrame({
        "test_id": test_ids,
        "converged": True,
        "degrees_of_freedom": 1,
        "n_subjects": 6,
        "statistic": np.array([1.0, 2.0, 4.0, 6.0, 8.0, 10.0]) ** 2,
        "mean_difference_norm": [0.2, 0.3, 0.5, 0.7, 0.8, 1.0],
        "p_value": [0.36, 0.16, 0.01, 0.002, 0.0005, 0.0001],
    })
    null = pd.DataFrame([
        {"test_id": test_id, "replicate": replicate, "p_value": p_value}
        for replicate, p_value in ((0, 0.2), (1, 0.7))
        for test_id in test_ids
    ])
    moderated, moderated_null = moderate_scalar_tests(table, null)
    assert moderated.variance_prior_df.notna().all()
    assert moderated.variance_prior.notna().all()
    assert np.isfinite(moderated.p_value).all()
    assert np.isfinite(moderated_null.p_value).all()
    assert not np.allclose(moderated.p_value, table.p_value)
    assert not np.allclose(moderated_null.p_value, null.p_value)


def test_complete_family_merge_counts_failed_hypotheses_in_bh(tmp_path, monkeypatch):
    from extra_scripts import merge_paired_path_test
    from tealeaf.sc.ds_benchmark import benjamini_hochberg

    shard = tmp_path / "shard_0"
    shard.mkdir()
    table = pd.DataFrame({"test_id": ["a", "b", "failed"], "converged": [True, True, False], "degrees_of_freedom": [1, 1, 1], "n_subjects": [6, 6, 0], "p_value": [.001, .2, 1.]})
    table["path_pseudocount"] = 32.
    table["path_pseudocount_scaling"] = "total"
    table.to_csv(shard / "paired_path.tsv", sep="\t", index=False)
    output = tmp_path / "merged"
    monkeypatch.setattr(merge_paired_path_test, "parse_args", lambda: SimpleNamespace(shards=[shard], output_dir=output, min_stratum_tests=100, calibration="native", moderate_variances=False, max_null_replicates=None, retain_failed_family=True))
    merge_paired_path_test.main()
    observed = pd.read_csv(output / "paired_path.tsv", sep="\t")
    assert np.allclose(observed.fdr, benjamini_hochberg(np.array([.001, .2, 1.])))
    assert observed.loc[observed.test_id.eq("failed"), "p_value"].iloc[0] == 1


def direct_leave_event_out_calibration(table, null):
    """Independent small-panel reference, remove own rows before sorting."""
    observed, reference = table.copy(), null.copy()
    observed["raw_p_value"], reference["raw_p_value"] = observed.p_value, reference.p_value
    reference["calibration_stratum"] = reference.test_id.map(observed.set_index("test_id").calibration_stratum)
    observed["p_value"], reference["p_value"] = np.nan, np.nan
    for output in (observed, reference):
        for index, row in output.iterrows():
            if pd.isna(row.calibration_stratum):
                continue
            pool = reference.loc[reference.calibration_stratum.eq(row.calibration_stratum)]
            if not len(pool):
                continue
            other = np.sort(pool.loc[~pool.test_id.eq(row.test_id), "raw_p_value"].to_numpy(float))
            output.at[index, "p_value"] = (1 + np.searchsorted(other, row.raw_p_value, side="right")) / (1 + len(other))
    return observed, reference


@pytest.mark.parametrize("seed", [0, 18, 932])
def test_batched_null_calibration_exactly_preserves_exclusion_ties_and_order(seed):
    rng = np.random.default_rng(seed)
    values = np.array([0., -0., np.nextafter(0., 1.), .001, .05, .3, 1., np.nan])
    table = pd.DataFrame({"test_id": [f"t{index}" for index in range(20)], "calibration_stratum": [f"s{index % 3}" for index in range(20)], "p_value": rng.choice(values, 20)})
    # Empty training stratum and a failed/no-own-null observation are retained.
    table.loc[19, "calibration_stratum"] = "empty"
    null = pd.DataFrame([{"test_id": f"t{index}", "replicate": replicate, "p_value": rng.choice(values)} for index in range(18) for replicate in range(int(rng.integers(1, 33)))])
    null = null.sample(frac=1., random_state=seed)
    null.index = np.arange(len(null)) * 3 + 7
    table.index = np.arange(len(table)) * 5 + 2
    expected = direct_leave_event_out_calibration(table, null)
    actual = empirical_null_calibration(table, null)
    for left, right in zip(actual, expected):
        pd.testing.assert_frame_equal(left, right, check_exact=True)


def test_single_event_null_pool_has_no_other_training_draws():
    table = pd.DataFrame({"test_id": ["only"], "calibration_stratum": ["single"], "p_value": [0.]})
    null = pd.DataFrame({"test_id": ["only"] * 4, "p_value": [0., .2, 1., np.nan]})
    observed, reference = empirical_null_calibration(table, null)
    np.testing.assert_array_equal(observed.p_value, [1.])
    np.testing.assert_array_equal(reference.p_value, np.ones(4))
