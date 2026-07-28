"""Numerical checks for the hierarchical multinomial model."""

import tempfile
from pathlib import Path
import unittest

import numpy as np
import scipy.sparse as sp

try:
    from tealeaf.sc import hierarchical

    hierarchical._torch()
except ImportError:  # pragma: no cover
    TORCH_AVAILABLE = False
else:
    TORCH_AVAILABLE = True


@unittest.skipUnless(TORCH_AVAILABLE, "PyTorch is unavailable")
class HierarchicalModelTest(unittest.TestCase):
    def setUp(self):
        design = sp.csr_matrix(
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.5, 0.5, 0.0],
                ],
                dtype=np.float32,
            )
        )
        counts = sp.csr_matrix(
            np.array(
                [[8, 2, 5, 1], [2, 8, 5, 1], [5, 5, 0, 2]],
                dtype=np.int64,
            )
        )
        gene_counts = sp.csr_matrix(
            np.array([[11, 5], [11, 5], [12, 0]], dtype=np.int64)
        )
        self.data = hierarchical.HierarchicalData(
            counts=(counts, counts.copy()),
            designs=(design, design.copy()),
            gene_counts=gene_counts,
            barcodes=np.array(["c1", "c2", "c3"]),
            transcripts=np.array(["t1", "t2", "t3"]),
            genes=np.array(["g1", "g2"]),
            transcript_gene=np.array([0, 0, 1]),
            primer_names=("poly", "hex"),
        )

    def test_within_gene_softmax_and_event_likelihood_match_dense(self):
        torch = hierarchical._torch()
        model = hierarchical.HierarchicalModel(
            self.data,
            rank=1,
            device="cpu",
            initial_state=np.array([[-1.0], [1.0], [0.0]], dtype=np.float32),
        )
        with torch.no_grad():
            model.gene_intercept.copy_(torch.tensor([0.2, -0.1]))
            model.gene_loadings.copy_(torch.tensor([[0.3], [-0.2]]))
            model.isoform_intercept.copy_(torch.tensor([0.1, -0.1, 0.0]))
            model.isoform_loadings.copy_(torch.tensor([[0.5], [-0.5], [0.0]]))
        observed = float(model.ec_loss(np.arange(3)).detach())

        state = model.state.detach().numpy()[:, 0]
        alpha = model.gene_intercept.detach().numpy()
        gene_loadings = model.gene_loadings.detach().numpy()[:, 0]
        beta = model.isoform_intercept.detach().numpy()
        iso_loadings = model.isoform_loadings.detach().numpy()[:, 0]
        expected = 0.0
        molecules = 0.0
        for cell in range(3):
            gene = np.exp(alpha + gene_loadings * state[cell])
            eta = beta + iso_loadings * state[cell]
            iso = np.array(
                [
                    np.exp(eta[0]) / np.exp(eta[:2]).sum(),
                    np.exp(eta[1]) / np.exp(eta[:2]).sum(),
                    1.0,
                ]
            )
            theta = gene[self.data.transcript_gene] * iso
            for counts, design in zip(self.data.counts, self.data.designs):
                prediction = design @ theta
                row = counts.getrow(cell).toarray().ravel()
                expected -= np.dot(row, np.log(prediction))
                expected += row.sum() * np.log(prediction.sum())
                molecules += row.sum()
        self.assertAlmostEqual(observed, expected / molecules, places=5)

    def test_fit_reduces_gene_loss_and_writes_factor_state(self):
        model = hierarchical.HierarchicalModel(
            self.data,
            rank=2,
            device="cpu",
            seed=2,
        )
        before = float(model.gene_loss(np.arange(3)).detach())
        history = model.fit_stage(
            "gene",
            epochs=15,
            batch_cells=3,
            learning_rate=0.05,
            seed=2,
        )
        after = float(model.gene_loss(np.arange(3)).detach())
        self.assertLess(after, before)
        self.assertEqual(len(history), 15)
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "fit_"
            model.write(prefix)
            with np.load(f"{prefix}glm_factors.npz") as saved:
                self.assertEqual(saved["right"].shape, (3, 2))
                self.assertEqual(saved["left"].shape, (3, 2))

    def test_isoform_stage_reduces_ec_loss(self):
        model = hierarchical.HierarchicalModel(
            self.data,
            rank=1,
            device="cpu",
            initial_state=np.array([[-1.0], [1.0], [0.0]], dtype=np.float32),
        )
        before = float(model.ec_loss(np.arange(3)).detach())
        model.fit_stage(
            "isoform",
            epochs=20,
            batch_cells=3,
            learning_rate=0.05,
            isoform_penalty=0.0,
            seed=3,
        )
        after = float(model.ec_loss(np.arange(3)).detach())
        self.assertLess(after, before)

    def test_sparse_gene_pca_is_finite(self):
        embedding, selected, components, diagnostics = (
            hierarchical.log_gene_randomized_pca(
                self.data.gene_counts,
                rank=1,
                n_hvg=2,
                power_iterations=1,
            )
        )
        self.assertEqual(embedding.shape, (3, 1))
        self.assertEqual(len(selected), 2)
        self.assertEqual(components.shape, (1, 2))
        self.assertTrue(np.isfinite(embedding).all())
        self.assertEqual(diagnostics["gene_pca_rank"], 1)
