import numpy as np

from extra_scripts.audit_path_measurement_covariance import summarize_bootstrap


def test_conditional_covariance_and_bias_have_distinct_denominators():
    proportions = np.asarray([[[.2, .8]], [[.4, .6]], [[.6, .4]]])
    covariance = np.asarray([[[.04, -.04], [-.04, .04]]])
    records = [{"subjects": ["subject"], "labels": [0], "proportions": value, "proportion_covariances": covariance, "scalar_proportion_covariances": covariance * 2, "effective_depths": [10.], "depths": [100.]} for value in proportions]
    output = summarize_bootstrap(records, {("subject", "0"): np.asarray([.3, .7])}, 4)
    assert len(output) == 1
    assert not output[0]["complete_bootstrap"]
    assert output[0]["bootstrap_converged"] == 3
    assert np.isclose(output[0]["empirical_trace"], .08)
    assert np.isclose(output[0]["fisher_variance_ratio"], 1)
    assert np.isclose(output[0]["scalar_variance_ratio"], 2)
    assert np.isclose(output[0]["bias_squared_norm"], .02)
    assert np.isclose(output[0]["mean_squared_error"], .08 * 2 / 3 + .02)
