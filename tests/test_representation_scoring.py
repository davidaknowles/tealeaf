"""Checks for reference-label scoring of compact cell factors."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

try:
    import sklearn  # noqa: F401
except ImportError:  # pragma: no cover - optional benchmark dependency
    SKLEARN_AVAILABLE = False
else:
    SKLEARN_AVAILABLE = True

from tealeaf.sc import representation_scoring


@unittest.skipUnless(SKLEARN_AVAILABLE, "scikit-learn is unavailable")
class RepresentationScoringTest(unittest.TestCase):
    def test_blocked_log_gene_reconstruction_matches_dense(self):
        left = np.array(
            [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]], dtype=np.float32
        )
        right = np.array(
            [[1.0, 0.2], [0.2, 1.0], [0.7, 0.3], [0.1, 0.8]],
            dtype=np.float32,
        )
        mapping = pd.DataFrame(
            {
                "transcript_id": ["t1", "t2", "t3"],
                "gene_id": ["g1", "g1", "g2"],
            }
        )
        gene_left, diagnostics = representation_scoring.aggregate_transcript_loadings(
            left, ["t1", "t2", "t3"], mapping, ["g1", "g2", "g3"]
        )
        expected_gene_left = np.array(
            [[1.5, 0.5], [0.0, 1.0], [0.0, 0.0]], dtype=np.float32
        )
        np.testing.assert_allclose(gene_left, expected_gene_left)
        self.assertEqual(diagnostics["n_mapped_transcripts"], 3)

        expression = representation_scoring.FactorizedGeneExpression(
            right, gene_left, target_sum=100.0, batch_cells=3, device="cpu"
        )
        blocked = np.vstack([block for _, block in expression.iter_blocks()])
        dense = np.maximum(right @ expected_gene_left.T, 0)
        dense = np.log1p(dense * (100.0 / dense.sum(axis=1))[:, None])
        np.testing.assert_allclose(blocked, dense, rtol=1e-6, atol=1e-6)
        means, variances = expression.mean_variance()
        np.testing.assert_allclose(means, dense.mean(axis=0), atol=1e-6)
        np.testing.assert_allclose(variances, dense.var(axis=0, ddof=1), atol=1e-6)

    def test_log_gene_incremental_pca_matches_dense_geometry(self):
        rng = np.random.default_rng(3)
        right = rng.uniform(0.1, 1.0, size=(20, 3)).astype(np.float32)
        gene_left = rng.uniform(0.1, 1.0, size=(6, 3)).astype(np.float32)
        embedding, genes, active, diagnostics = (
            representation_scoring.log_gene_pca_embedding(
                right,
                gene_left,
                np.array([f"g{i}" for i in range(6)]),
                target_sum=100.0,
                n_hvg=4,
                n_components=3,
                batch_cells=32,
                device="cpu",
            )
        )
        dense = right @ gene_left.T
        dense = np.log1p(dense * (100.0 / dense.sum(axis=1))[:, None])
        selected = np.array([int(gene[1:]) for gene in genes])
        from sklearn.decomposition import PCA

        expected = PCA(n_components=3).fit_transform(dense[:, selected])
        np.testing.assert_allclose(
            squareform(pdist(embedding)),
            squareform(pdist(expected)),
            rtol=1e-4,
            atol=1e-5,
        )
        self.assertTrue(active.all())
        self.assertEqual(diagnostics["n_hvg"], 4)
        self.assertEqual(diagnostics["representation"], "log1p_gene_pca")
        self.assertGreaterEqual(diagnostics["pca_explained_variance"], 0.0)
        self.assertLessEqual(diagnostics["pca_explained_variance"], 1.0)

    def test_log_gene_pca_rejects_library_size_only_factor(self):
        right = np.arange(1, 21, dtype=np.float32)[:, None]
        gene_left = np.array([[1.0], [2.0], [4.0]], dtype=np.float32)
        with self.assertRaisesRegex(ValueError, "representation is collapsed"):
            representation_scoring.log_gene_pca_embedding(
                right,
                gene_left,
                np.array(["g1", "g2", "g3"]),
                target_sum=100.0,
                n_hvg=3,
                n_components=2,
                batch_cells=8,
                device="cpu",
            )

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
