import numpy as np
import pandas as pd

from extra_scripts.merge_paired_path_test import moderate_scalar_tests


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
