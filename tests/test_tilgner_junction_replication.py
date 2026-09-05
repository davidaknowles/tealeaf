import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.assess_tilgner_junction_replication import junction_from_majiq, mean_composition, mean_composition_batch, paired_logratio_two_category_batch, replace_majiq_results
from tealeaf.sc.differential import paired_logratio_test


def test_mean_composition_weights_pseudobulks_equally():
    counts = np.array([[9, 1], [1, 1]])
    assert np.allclose(mean_composition(counts), [0.7, 0.3])


def test_batched_mean_composition_matches_scalar_version():
    counts = sparse.csr_matrix([[9, 1, 4, 6], [1, 1, 3, 7]])
    observed = mean_composition_batch(counts, np.arange(2), [[0, 1], [2, 3]])
    assert np.allclose(observed[0], mean_composition(counts[:, [0, 1]]))
    assert np.allclose(observed[1], mean_composition(counts[:, [2, 3]]))


def test_majiq_boundaries_become_intron_coordinates():
    row = type("Row", (), {"seqid": "chr1", "start": 100, "end": 200, "strand": "+"})()
    assert junction_from_majiq(row) == ("chr1", 101, 199, "+")


def test_batched_two_category_clr_matches_scalar_test():
    counts = np.array([[8, 2, 3, 7], [5, 5, 4, 6], [9, 1, 2, 8], [6, 4, 5, 5], [7, 3, 1, 9], [4, 6, 6, 4], [10, 0, 4, 6], [7, 3, 7, 3]])
    labels = np.tile([0, 1], 4)
    clusters = np.repeat(np.arange(4), 2)
    batched = paired_logratio_two_category_batch(sparse.csr_matrix(counts), np.arange(8), [[0, 1], [2, 3]], minimum_subjects=3)
    for indices, (p_value, n_subjects) in zip(([0, 1], [2, 3]), batched):
        scalar = paired_logratio_test(counts[:, indices], labels, clusters)
        assert np.isclose(p_value, scalar["p_value"])
        assert n_subjects == scalar["n_subjects"]


def test_replace_majiq_results(tmp_path):
    columns = ["method", "contrast_id", "effect", "level_a", "level_b", "significant"]
    original = pd.DataFrame([["LeafCutter", "cell_type__A__B", "cell_type", "A", "B", True], ["MAJIQ Heterogen", "cell_type__A__B", "cell_type", "A", "B", False]], columns=columns)
    replacement = pd.DataFrame([["MAJIQ Heterogen", "cell_type__A__B", "cell_type", "A", "B", True]], columns=columns)
    replacement.to_csv(tmp_path / "cell_type__A__B.tsv", sep="\t", index=False)
    observed = replace_majiq_results(original, tmp_path, {"A", "B"})
    assert observed["method"].tolist() == ["LeafCutter", "MAJIQ Heterogen"]
    assert observed.loc[observed["method"] == "MAJIQ Heterogen", "significant"].item()
