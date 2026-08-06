import numpy as np
import pandas as pd

from extra_scripts.select_ec_glmm_hybrid import select_hybrid


def test_hybrid_selects_noise_only_when_null_bic_improves():
    base = pd.DataFrame({
        "test_id": ["a", "b"],
        "method": "plain",
        "n_samples": [20, 20],
        "null_objective": [100.0, 100.0],
        "null_converged": True,
        "alternative_converged": True,
        "p_value": [0.01, 0.02],
        "fdr": np.nan,
    })
    noise = base.copy()
    noise["method"] = "noise"
    noise["null_objective"] = [95.0, 99.9]
    noise["p_value"] = [0.03, 0.04]
    result = select_hybrid(base, noise).set_index("test_id")
    assert result.loc["a", "selected_observation_model"] == "logistic_normal_noise"
    assert result.loc["a", "p_value"] == 0.03
    assert result.loc["b", "selected_observation_model"] == "multinomial"
    assert result.loc["b", "p_value"] == 0.02
    assert result.method.eq("laplace_multinomial_hybrid_bic").all()
