"""Checks for reference-label scoring of compact cell factors."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import sklearn  # noqa: F401
except ImportError:  # pragma: no cover - optional benchmark dependency
    SKLEARN_AVAILABLE = False
else:
    SKLEARN_AVAILABLE = True

from tealeaf.sc import representation_scoring


@unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn is unavailable")
class RepresentationScoringTest(unittest.TestCase):
    def test_grouped_scores_recover_separated_labels(self):
        rng = np.random.default_rng(0)
        labels = np.tile(np.repeat(["a", "b", "c"], 5), 4)
        groups = np.repeat(["m1", "m2", "m3", "m4"], 15)
        centers = {label: np.eye(3)[index] for index, label in enumerate(["a", "b", "c"])}
        representation = np.vstack([centers[label] for label in labels])
        representation += rng.normal(0, 0.01, representation.shape)
        report, folds = representation_scoring.score_representation(
            representation, labels, groups, n_splits=4
        )
        self.assertEqual(len(folds), 4)
        self.assertGreater(report["macro_f1_mean"], 0.95)
        self.assertGreater(report["adjusted_rand_index"], 0.95)
        self.assertGreater(report["reference_label_silhouette"], 0.9)
        self.assertGreater(report["kmeans_cluster_silhouette"], 0.9)
        self.assertEqual(report["pca_components"], 3)

    def test_load_and_align_factors(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "fit_"
            np.savetxt(f"{prefix}glm_rows.txt", ["c1", "c2", "c3"], fmt="%s")
            np.savez_compressed(
                f"{prefix}glm_factors.npz",
                left=np.ones((2, 2)),
                right=np.arange(6).reshape(3, 2),
            )
            cell_ids, factors = representation_scoring.load_factor_representation(prefix)
            labels = pd.DataFrame(
                {"cell_id": ["c3", "c1"], "label": ["b", "a"]}
            )
            aligned, aligned_labels, _, aligned_ids = (
                representation_scoring.align_reference_labels(
                    cell_ids, factors, labels
                )
            )
            np.testing.assert_array_equal(aligned_ids, ["c3", "c1"])
            np.testing.assert_array_equal(aligned_labels, ["b", "a"])
            np.testing.assert_array_equal(aligned, [[4, 5], [0, 1]])


if __name__ == "__main__":
    unittest.main()
