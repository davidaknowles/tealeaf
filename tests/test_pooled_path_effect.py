import numpy as np
import pytest

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.ec_block_glmm import pooled_path_effect


def test_pooled_effect_has_b_minus_a_orientation_and_preserves_inputs():
    first = np.array([[90., 10.], [10., 90.]])
    second = np.array([[2., 18.], [18., 2.]])
    data = ECGLMMData((first, second), (np.eye(2), np.eye(2)), np.ones((2, 1)), np.array(["s", "s"]))
    baseline = np.array([.5, .5])
    pooled = pooled_path_effect(data, [0, 1], [0, 1], baseline=baseline)
    balanced = pooled_path_effect(data, [0, 1], [0, 1], baseline=baseline, balance_primers=True)
    ran = pooled_path_effect(data, [0, 1], [0, 1], baseline=baseline, primer=0)
    oligo = pooled_path_effect(data, [0, 1], [0, 1], baseline=baseline, primer=1)
    assert pooled["converged"] and pooled["difference"][1] > 0
    assert np.linalg.norm(balanced["difference"]) < 1e-6
    assert ran["difference"][1] > 0 and oligo["difference"][1] < 0
    assert np.array_equal(first, [[90., 10.], [10., 90.]])
    assert np.array_equal(second, [[2., 18.], [18., 2.]])
    assert np.array_equal(baseline, [.5, .5])


def test_pooled_effect_rejects_unaligned_labels_or_invalid_primer():
    data = ECGLMMData((np.ones((2, 2)),), (np.eye(2),), np.ones((2, 1)), np.array(["s", "s"]))
    with pytest.raises(ValueError, match="aligned"):
        pooled_path_effect(data, [0, 1], [0])
    with pytest.raises(ValueError, match="primer"):
        pooled_path_effect(data, [0, 1], [0, 1], primer=2)
