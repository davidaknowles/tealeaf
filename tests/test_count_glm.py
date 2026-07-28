"""Numerical checks for the count-aware flat GLM."""

import unittest

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import count_glm
from tealeaf.sc.hierarchical import HierarchicalData


class CountAwareGLMTest(unittest.TestCase):
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
        self.data = HierarchicalData(
            counts=(counts, counts.copy()),
            designs=(design, design.copy()),
            gene_counts=sp.csr_matrix(
                np.array([[11, 5], [11, 5], [12, 0]], dtype=np.int64)
            ),
            barcodes=np.array(["c1", "c2", "c3"]),
            transcripts=np.array(["t1", "t2", "t3"]),
            genes=np.array(["g1", "g2"]),
            transcript_gene=np.array([0, 0, 1]),
            primer_names=("poly", "hex"),
        )
        self.model = count_glm.CountAwareGLM(
            self.data,
            rank=1,
            selected_genes=np.array([0, 1]),
            concentration=5.0,
            device="cpu",
            initial_state=np.array([[-1.0], [1.0], [0.0]], dtype=np.float32),
        )

    def test_nb_loss_matches_scipy_distribution(self):
        from scipy.stats import nbinom

        torch = count_glm._torch()
        cells = np.arange(3)
        observed_loss = float(self.model.gene_loss(cells).detach())
        with torch.no_grad():
            theta = self.model._abundance(torch.as_tensor(cells, dtype=torch.int64))
            mass = torch.segment_reduce(
                theta,
                reduce="sum",
                lengths=self.model.gene_lengths,
                axis=0,
            ).numpy()
        observed = self.data.gene_counts.toarray().T
        totals = observed.sum(axis=0)
        expected = mass / mass.sum(axis=0) * totals
        concentration = self.model.concentration
        probability = concentration / (concentration + expected)
        log_probability = nbinom.logpmf(observed, concentration, probability)
        self.assertAlmostEqual(
            observed_loss,
            -float(log_probability.sum()) / totals.sum(),
            places=5,
        )

    def test_gene_fit_decreases_nb_loss(self):
        cells = np.arange(3)
        before = float(self.model.gene_loss(cells).detach())
        history = self.model.fit_stage(
            "gene",
            epochs=20,
            batch_cells=3,
            learning_rate=0.05,
            seed=4,
        )
        after = float(self.model.gene_loss(cells).detach())
        self.assertLess(after, before)
        self.assertEqual(len(history), 20)

    def test_poisson_loss_matches_scipy_distribution(self):
        from scipy.stats import poisson

        model = count_glm.CountAwareGLM(
            self.data,
            rank=1,
            selected_genes=np.array([0, 1]),
            gene_family="poisson",
            device="cpu",
            initial_state=np.array([[-1.0], [1.0], [0.0]], dtype=np.float32),
        )
        torch = count_glm._torch()
        cells = np.arange(3)
        observed_loss = float(model.gene_loss(cells).detach())
        with torch.no_grad():
            theta = model._abundance(torch.as_tensor(cells, dtype=torch.int64))
            mass = torch.segment_reduce(
                theta,
                reduce="sum",
                lengths=model.gene_lengths,
                axis=0,
            ).numpy()
        observed = self.data.gene_counts.toarray().T
        totals = observed.sum(axis=0)
        expected = mass / mass.sum(axis=0) * totals
        self.assertAlmostEqual(
            observed_loss,
            -float(poisson.logpmf(observed, expected).sum()) / totals.sum(),
            places=5,
        )

    def test_standardized_log_gene_fit_decreases_loss(self):
        model = count_glm.CountAwareGLM(
            self.data,
            rank=1,
            selected_genes=np.array([0, 1]),
            gene_family="standardized_log",
            device="cpu",
            initial_state=np.array([[-1.0], [1.0], [0.0]], dtype=np.float32),
        )
        cells = np.arange(3)
        before = float(model.gene_loss(cells).detach())
        model.fit_stage(
            "gene",
            epochs=20,
            batch_cells=3,
            learning_rate=0.05,
            seed=5,
        )
        after = float(model.gene_loss(cells).detach())
        self.assertLess(after, before)
