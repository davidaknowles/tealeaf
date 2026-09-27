import numpy as np

from extra_scripts.run_suppa2_tealeaf_hybrid import (
    collapse_event_nuisance,
    event_path_index,
    partition_event_tests,
)
from tealeaf.sc.ec_glmm import ECGLMMData


def test_event_path_index_preserves_non_event_isoforms_as_nuisance():
    features = ["txA.1", "txB.2", "txC.1", "txD.4"]
    result = event_path_index(
        np.array([0, 1, 2, 3]),
        features,
        "txA.1,txC.9",
        "txB.2",
    )
    np.testing.assert_array_equal(result, [0, 1, 0, -1])


def test_event_path_index_rejects_incomplete_or_ambiguous_event():
    features = ["txA.1", "txB.1"]
    assert event_path_index([0, 1], features, "txA.1", "txMissing.1") is None
    assert event_path_index([0, 1], features, "txA.1", "txA.2") is None


def test_partition_event_tests_keeps_context_together_and_balances():
    def candidate(gene, rows):
        return ("", "", "", gene, [], [], [], np.array(rows), None, ("A", "B"))

    tests = [
        (candidate(0, [0, 1]), "e1", None),
        (candidate(0, [0, 1]), "e2", None),
        (candidate(1, [2, 3]), "e3", None),
        (candidate(2, [4, 5]), "e4", None),
    ]
    shards = partition_event_tests(tests, 2)
    assert sorted(map(len, shards)) == [2, 2]
    contexts = [
        {(test[0][3], tuple(test[0][7])) for test in shard}
        for shard in shards
    ]
    assert not contexts[0] & contexts[1]


def test_collapse_event_nuisance_preserves_fixed_mixture_and_event_paths():
    mapping = np.array(
        [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 0, 1]], dtype=float
    )
    data = ECGLMMData(
        counts=(np.zeros((2, 3)),),
        compatibility=(mapping,),
        design=np.ones((2, 1)),
        clusters=np.array(["s1", "s2"]),
    )
    collapsed, path_index, baseline = collapse_event_nuisance(
        data, np.array([0, 1, -1, -1]), np.array([0.1, 0.2, 0.3, 0.4])
    )
    np.testing.assert_array_equal(path_index, [0, 1, -1])
    np.testing.assert_allclose(baseline, [0.1, 0.2, 0.7])
    expected_nuisance = (0.3 * mapping[:, 2] + 0.4 * mapping[:, 3]) / 0.7
    np.testing.assert_allclose(collapsed.compatibility[0][:, 0], mapping[:, 0])
    np.testing.assert_allclose(collapsed.compatibility[0][:, 1], mapping[:, 1])
    np.testing.assert_allclose(collapsed.compatibility[0][:, 2], expected_nuisance)
