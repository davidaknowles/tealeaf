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


def test_mixed_event_record_profiles_transcripts_and_exports_complete_psi():
    from types import SimpleNamespace
    from extra_scripts.run_suppa2_tealeaf_hybrid import mixed_event_record
    subjects = np.repeat(np.arange(6), 2)
    labels = np.tile([0, 1], 6)
    counts = np.tile([[30., 20., 40., 10.], [50., 20., 20., 10.]], (6, 1))
    base = ECGLMMData((counts,), (np.eye(4),), np.ones((12, 1)), subjects)
    event = SimpleNamespace(event_id="gene;SE:example", feature_id="SUPPA2:gene;SE:example", event_type="SE")
    args = SimpleNamespace(max_iter=300, null_multistart=False, report_pseudocount=1., null_replicates=2, seed=20260927, export_path_usage=True)
    archive = ([], [])
    record, null, usage = mixed_event_record(base, np.array([0, 0, 1, -1]), labels, subjects, np.array([.4, .2, .3, .1]), event, "gene", ("A", "B"), counts.sum(axis=1), 4, args, score_archive=archive)
    assert record["converged"] and record["complete_reporting_fits"]
    assert record["n_isoforms"] == record["n_source_isoforms"] == 4
    assert record["n_expected_subjects"] == record["n_fitted_subjects"] == 6
    assert record["path_pseudocount"] == 0 and record["report_pseudocount"] == 1
    assert record["effect_size"] > 0 and record["test_ilr_effect_size"] > 0
    assert len(null) == 2 and len(usage) == 12
    assert {row["replicate"] for row in null} == {0, 1}
    assert len(archive[0]) == 1 and len(archive[1]) == 6
    assert archive[0][0]["n_expected_subjects"] == 6
    assert all(row["reference_information"] >= row["information"] for row in archive[1])
    from tealeaf.sc.path_score_mixed import shared_path_score_components
    components = shared_path_score_components(base, [0, 0, 1, -1], labels, subjects, baseline=np.array([.4, .2, .3, .1]), reporting_concentration=1.)
    reused, reused_null, reused_usage = mixed_event_record(base, np.array([0, 0, 1, -1]), labels, subjects, np.array([.4, .2, .3, .1]), event, "gene", ("A", "B"), counts.sum(axis=1), 4, args, components=components)
    assert reused == record and reused_null == null and reused_usage == usage
