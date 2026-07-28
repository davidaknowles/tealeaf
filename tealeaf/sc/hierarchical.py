"""Hierarchical gene/isoform multinomial model for paired single-cell EC counts."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import time

import numpy as np
import scipy.sparse as sp


def _torch():
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("The hierarchical model requires PyTorch") from exc
    return torch


@dataclass
class HierarchicalData:
    """Raw paired counts and observation matrices in one transcript ordering."""

    counts: tuple[sp.csr_matrix, sp.csr_matrix]
    designs: tuple[sp.csr_matrix, sp.csr_matrix]
    gene_counts: sp.csr_matrix
    barcodes: np.ndarray
    transcripts: np.ndarray
    genes: np.ndarray
    transcript_gene: np.ndarray
    primer_names: tuple[str, str]


def prepare_hierarchical_data(prepared, selected, transcript_to_gene):
    """Convert paired GLM preparation into raw-count hierarchical inputs."""
    from tealeaf.sc import glm_cv

    if prepared.cv_raw_counts is None:
        raise ValueError("hierarchical fitting requires paired raw counts")
    selected = np.asarray(selected, dtype=np.int64)
    raw = prepared.cv_raw_counts[selected].tocsr()
    if raw.shape[1] % 2:
        raise ValueError("paired raw counts require two equal-width blocks")
    n_ec = raw.shape[1] // 2
    if prepared.compatibility.shape[0] != 2 * n_ec:
        raise ValueError(
            "hierarchical fitting does not accept auxiliary response blocks"
        )
    designs = tuple(
        (2.0 * prepared.compatibility[start : start + n_ec]).tocsr()
        for start in (0, n_ec)
    )
    counts = (raw[:, :n_ec].tocsr(), raw[:, n_ec:].tocsr())
    transcript_gene_matrix, genes = glm_cv._transcript_gene_assignment(
        prepared.features, transcript_to_gene
    )
    mapped = np.asarray(transcript_gene_matrix.getnnz(axis=1)).ravel() == 1
    if not mapped.all():
        raise ValueError(
            f"{int((~mapped).sum())} retained transcripts lack a unique gene"
        )
    transcript_gene = np.asarray(transcript_gene_matrix.indices, dtype=np.int64)
    support = designs[0].copy()
    support.data.fill(1.0)
    other = designs[1].copy()
    other.data.fill(1.0)
    support = (support + other).tocsr()
    support.data.fill(1.0)
    ec_gene = glm_cv._unambiguous_ec_gene_assignment(support, transcript_gene_matrix)
    gene_counts = (counts[0] @ ec_gene + counts[1] @ ec_gene).tocsr()
    gene_counts.sum_duplicates()
    metadata = prepared.metadata or {}
    primer_names = tuple(metadata.get("primer_names", ("polydt", "ranhex")))
    return HierarchicalData(
        counts=counts,
        designs=designs,
        gene_counts=gene_counts,
        barcodes=np.asarray(prepared.barcodes)[selected],
        transcripts=np.asarray(prepared.features, dtype=str),
        genes=np.asarray(genes, dtype=str),
        transcript_gene=transcript_gene,
        primer_names=primer_names,
    )


def log_gene_randomized_pca(
    counts,
    rank,
    *,
    n_hvg=2_000,
    target_sum=10_000.0,
    oversamples=10,
    power_iterations=2,
    seed=0,
):
    """Return a centered/scaled randomized PCA without densifying gene counts."""
    counts = counts.astype(np.float32).tocsr(copy=True)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    active = totals > 0
    if active.sum() <= int(rank):
        raise ValueError("not enough nonempty cells for gene PCA")
    scales = np.zeros_like(totals, dtype=np.float32)
    scales[active] = float(target_sum) / totals[active]
    counts.data *= np.repeat(scales, np.diff(counts.indptr))
    np.log1p(counts.data, out=counts.data)
    means = np.asarray(counts.mean(axis=0)).ravel().astype(np.float64)
    second = np.asarray(counts.power(2).mean(axis=0)).ravel()
    variances = np.maximum(second - np.square(means), 0.0)
    n_hvg = min(max(int(rank), int(n_hvg)), counts.shape[1])
    selected = np.argsort(variances)[-n_hvg:]
    matrix = counts[:, selected].tocsr()
    means = means[selected]
    ddof = matrix.shape[0] / max(matrix.shape[0] - 1, 1)
    feature_scales = np.sqrt(np.maximum(variances[selected] * ddof, 0.0))
    feature_scales[~np.isfinite(feature_scales) | (feature_scales == 0)] = 1.0
    width = min(int(rank) + int(oversamples), matrix.shape[1])
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((matrix.shape[1], width))

    def right_multiply(vectors):
        scaled = vectors / feature_scales[:, None]
        return np.asarray(matrix @ scaled) - np.outer(
            np.ones(matrix.shape[0]), means @ scaled
        )

    def left_multiply(vectors):
        cross = np.asarray(matrix.T @ vectors)
        return (cross - np.outer(means, vectors.sum(axis=0))) / feature_scales[:, None]

    q, _ = np.linalg.qr(right_multiply(omega), mode="reduced")
    for _ in range(int(power_iterations)):
        q, _ = np.linalg.qr(right_multiply(left_multiply(q)), mode="reduced")
    small = left_multiply(q).T
    left, singular_values, right_t = np.linalg.svd(small, full_matrices=False)
    rank = min(int(rank), len(singular_values))
    embedding = (q @ left[:, :rank] * singular_values[:rank]).astype(np.float32)
    embedding[~active] = 0
    diagnostics = {
        "gene_pca_hvg": int(n_hvg),
        "gene_pca_rank": int(rank),
        "gene_pca_active_cells": int(active.sum()),
        "gene_pca_explained_variance": float(
            np.square(singular_values[:rank]).sum()
            / max(np.square(singular_values).sum(), np.finfo(float).tiny)
        ),
    }
    return embedding, selected, right_t[:rank].astype(np.float32), diagnostics


def _event_edges(counts, design):
    """Expand observed cell/EC events to their sparse design edges."""
    coo = counts.tocoo(copy=False)
    event_cells = np.asarray(coo.row, dtype=np.int64)
    event_ecs = np.asarray(coo.col, dtype=np.int64)
    event_counts = np.asarray(coo.data, dtype=np.float32)
    lengths = np.diff(design.indptr)[event_ecs]
    keep = lengths > 0
    event_cells = event_cells[keep]
    event_ecs = event_ecs[keep]
    event_counts = event_counts[keep]
    lengths = lengths[keep]
    if not len(lengths):
        empty = np.empty(0, dtype=np.int64)
        return (
            event_cells,
            event_counts,
            empty,
            empty,
            empty,
            np.empty(0, dtype=np.float32),
        )
    offsets = np.cumsum(lengths, dtype=np.int64) - lengths
    starts = design.indptr[event_ecs]
    positions = np.arange(lengths.sum(), dtype=np.int64)
    positions += np.repeat(starts - offsets, lengths)
    edge_events = np.repeat(np.arange(len(lengths), dtype=np.int64), lengths)
    edge_cells = np.repeat(event_cells, lengths)
    return (
        event_cells,
        event_counts,
        edge_events,
        edge_cells,
        np.asarray(design.indices[positions], dtype=np.int64),
        np.asarray(design.data[positions], dtype=np.float32),
    )


class HierarchicalModel:
    """One latent cell state with gene and within-gene isoform decoder heads."""

    def __init__(
        self,
        data,
        rank,
        *,
        device="auto",
        initial_state=None,
        seed=0,
    ):
        torch = _torch()
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.torch = torch
        self.device = torch.device(device)
        self.data = data
        self.rank = int(rank)
        self.n_cells = len(data.barcodes)
        self.n_genes = len(data.genes)
        self.n_transcripts = len(data.transcripts)
        order = np.argsort(data.transcript_gene, kind="stable")
        self.transcript_order = torch.as_tensor(
            order, dtype=torch.int64, device=self.device
        )
        inverse = np.empty_like(order)
        inverse[order] = np.arange(len(order))
        self.original_to_sorted = torch.as_tensor(
            inverse, dtype=torch.int64, device=self.device
        )
        sorted_gene = data.transcript_gene[order]
        self.sorted_gene = torch.as_tensor(
            sorted_gene, dtype=torch.int64, device=self.device
        )
        self.gene_lengths = torch.as_tensor(
            np.bincount(sorted_gene, minlength=self.n_genes),
            dtype=torch.int64,
            device=self.device,
        )
        rng = np.random.default_rng(seed)
        if initial_state is None:
            initial_state = rng.normal(0, 0.01, (self.n_cells, self.rank)).astype(
                np.float32
            )
        initial_state = np.asarray(initial_state, dtype=np.float32)
        if initial_state.shape != (self.n_cells, self.rank):
            raise ValueError("initial state has the wrong shape")
        initial_state -= initial_state.mean(axis=0, keepdims=True)
        state_sd = initial_state.std(axis=0, keepdims=True)
        state_sd[state_sd == 0] = 1
        initial_state /= state_sd
        gene_totals = np.asarray(data.gene_counts.sum(axis=0)).ravel()
        alpha = np.log(gene_totals + 0.5)
        alpha -= alpha.mean()
        self.state = torch.nn.Parameter(
            torch.as_tensor(initial_state, device=self.device)
        )
        self.gene_intercept = torch.nn.Parameter(
            torch.as_tensor(alpha, dtype=torch.float32, device=self.device)
        )
        self.gene_loadings = torch.nn.Parameter(
            torch.zeros(
                (self.n_genes, self.rank),
                dtype=torch.float32,
                device=self.device,
            )
        )
        self.isoform_intercept = torch.nn.Parameter(
            torch.zeros(self.n_transcripts, dtype=torch.float32, device=self.device)
        )
        self.isoform_loadings = torch.nn.Parameter(
            torch.zeros(
                (self.n_transcripts, self.rank),
                dtype=torch.float32,
                device=self.device,
            )
        )
        self.column_sums = tuple(
            torch.as_tensor(
                np.asarray(design.sum(axis=0)).ravel()[order],
                dtype=torch.float32,
                device=self.device,
            )
            for design in data.designs
        )

    def parameters(self, stage="joint"):
        if stage == "gene":
            return [self.state, self.gene_intercept, self.gene_loadings]
        if stage == "isoform":
            return [self.isoform_intercept, self.isoform_loadings]
        if stage == "joint":
            return [
                self.state,
                self.gene_intercept,
                self.gene_loadings,
                self.isoform_intercept,
                self.isoform_loadings,
            ]
        raise ValueError("stage must be gene, isoform, or joint")

    def _abundance(self, cell_indices):
        torch = self.torch
        state = self.state[cell_indices]
        gene_logits = (
            self.gene_intercept[:, None] + self.gene_loadings @ state.T
        ).clamp(-20, 20)
        gene_abundance = torch.exp(gene_logits)
        eta = self.isoform_intercept[:, None] + self.isoform_loadings @ state.T
        eta = eta.index_select(0, self.transcript_order)
        maxima = torch.segment_reduce(
            eta, reduce="max", lengths=self.gene_lengths, axis=0
        )
        shifted = eta - torch.repeat_interleave(maxima, self.gene_lengths, dim=0)
        unnormalized = torch.exp(shifted)
        sums = torch.segment_reduce(
            unnormalized, reduce="sum", lengths=self.gene_lengths, axis=0
        )
        proportions = unnormalized / torch.repeat_interleave(
            sums.clamp_min(1e-30), self.gene_lengths, dim=0
        )
        theta = proportions * gene_abundance.index_select(0, self.sorted_gene)
        return theta, gene_logits

    def gene_loss(self, cell_indices):
        torch = self.torch
        cell_indices = np.asarray(cell_indices, dtype=np.int64)
        block = self.data.gene_counts[cell_indices].tocoo(copy=False)
        indices = torch.as_tensor(cell_indices, dtype=torch.int64, device=self.device)
        state = self.state[indices]
        logits = (self.gene_intercept[:, None] + self.gene_loadings @ state.T).clamp(
            -20, 20
        )
        totals = np.asarray(self.data.gene_counts[cell_indices].sum(axis=1)).ravel()
        total = max(float(totals.sum()), 1.0)
        loss = (
            torch.as_tensor(totals, dtype=torch.float32, device=self.device)
            * torch.logsumexp(logits, dim=0)
        ).sum()
        if block.nnz:
            rows = torch.as_tensor(block.row, dtype=torch.int64, device=self.device)
            cols = torch.as_tensor(block.col, dtype=torch.int64, device=self.device)
            values = torch.as_tensor(
                block.data, dtype=torch.float32, device=self.device
            )
            loss -= (values * logits[cols, rows]).sum()
        return loss / total

    def ec_loss(self, cell_indices, isoform_penalty=0.0):
        torch = self.torch
        cell_indices = np.asarray(cell_indices, dtype=np.int64)
        tensor_cells = torch.as_tensor(
            cell_indices, dtype=torch.int64, device=self.device
        )
        theta, _ = self._abundance(tensor_cells)
        total_loss = torch.zeros((), device=self.device)
        total_molecules = 0.0
        for counts, design, column_sums in zip(
            self.data.counts, self.data.designs, self.column_sums
        ):
            block = counts[cell_indices].tocsr()
            totals = np.asarray(block.sum(axis=1)).ravel().astype(np.float32)
            total_molecules += float(totals.sum())
            normalizer = (theta * column_sums[:, None]).sum(dim=0).clamp_min(1e-30)
            (
                event_cells,
                event_counts,
                edge_events,
                edge_cells,
                edge_transcripts,
                edge_weights,
            ) = _event_edges(block, design)
            if len(event_counts):
                edge_events_t = torch.as_tensor(
                    edge_events, dtype=torch.int64, device=self.device
                )
                edge_cells_t = torch.as_tensor(
                    edge_cells, dtype=torch.int64, device=self.device
                )
                sorted_transcripts = self.original_to_sorted.index_select(
                    0,
                    torch.as_tensor(
                        edge_transcripts,
                        dtype=torch.int64,
                        device=self.device,
                    ),
                )
                weights = torch.as_tensor(
                    edge_weights, dtype=torch.float32, device=self.device
                )
                predictions = torch.zeros(
                    len(event_counts), dtype=torch.float32, device=self.device
                )
                predictions.scatter_add_(
                    0,
                    edge_events_t,
                    weights * theta[sorted_transcripts, edge_cells_t],
                )
                event_cells_t = torch.as_tensor(
                    event_cells, dtype=torch.int64, device=self.device
                )
                values = torch.as_tensor(
                    event_counts, dtype=torch.float32, device=self.device
                )
                total_loss += (
                    values
                    * (
                        torch.log(normalizer[event_cells_t])
                        - torch.log(predictions.clamp_min(1e-30))
                    )
                ).sum()
        loss = total_loss / max(total_molecules, 1.0)
        if isoform_penalty:
            loss = loss + 0.5 * float(isoform_penalty) * (
                self.isoform_loadings.square().mean()
            )
        return loss

    def center_identifiable_parameters(self):
        """Apply likelihood-preserving intercept and within-gene centering."""
        torch = self.torch
        with torch.no_grad():
            self.gene_intercept.sub_(self.gene_intercept.mean())
            beta = self.isoform_intercept.index_select(0, self.transcript_order)
            beta_means = torch.segment_reduce(
                beta, reduce="mean", lengths=self.gene_lengths, axis=0
            )
            beta -= torch.repeat_interleave(beta_means, self.gene_lengths, dim=0)
            self.isoform_intercept.index_copy_(0, self.transcript_order, beta)
            loadings = self.isoform_loadings.index_select(0, self.transcript_order)
            loading_means = torch.segment_reduce(
                loadings, reduce="mean", lengths=self.gene_lengths, axis=0
            )
            loadings -= torch.repeat_interleave(loading_means, self.gene_lengths, dim=0)
            self.isoform_loadings.index_copy_(0, self.transcript_order, loadings)

    def fit_stage(
        self,
        stage,
        *,
        epochs,
        batch_cells,
        learning_rate,
        isoform_penalty=0.0,
        joint_gene_weight=1.0,
        seed=0,
        progress_callback=None,
    ):
        """Fit one stage with shuffled cell minibatches."""
        torch = self.torch
        optimizer = torch.optim.Adam(self.parameters(stage), lr=float(learning_rate))
        rng = np.random.default_rng(seed)
        history = []
        for epoch in range(int(epochs)):
            started = time.perf_counter()
            order = rng.permutation(self.n_cells)
            weighted_loss = 0.0
            weighted_gene_loss = 0.0
            weighted_ec_loss = 0.0
            seen = 0
            for start in range(0, self.n_cells, int(batch_cells)):
                batch = order[start : start + int(batch_cells)]
                optimizer.zero_grad(set_to_none=True)
                gene_component = None
                ec_component = None
                if stage == "gene":
                    gene_component = self.gene_loss(batch)
                    loss = gene_component
                else:
                    ec_component = self.ec_loss(batch, isoform_penalty)
                    loss = ec_component
                    if stage == "joint" and joint_gene_weight:
                        gene_component = self.gene_loss(batch)
                        loss = loss + float(joint_gene_weight) * gene_component
                if not bool(torch.isfinite(loss).item()):
                    raise FloatingPointError(
                        f"non-finite {stage} loss at epoch {epoch + 1}"
                    )
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.parameters(stage), max_norm=10.0)
                optimizer.step()
                weighted_loss += float(loss.detach().cpu()) * len(batch)
                seen += len(batch)
                if gene_component is not None:
                    weighted_gene_loss += float(gene_component.detach().cpu()) * len(
                        batch
                    )
                if ec_component is not None:
                    weighted_ec_loss += float(ec_component.detach().cpu()) * len(batch)
            self.center_identifiable_parameters()
            row = {
                "stage": stage,
                "epoch": epoch + 1,
                "loss": weighted_loss / max(seen, 1),
                "seconds": time.perf_counter() - started,
            }
            if weighted_gene_loss:
                row["gene_loss"] = weighted_gene_loss / max(seen, 1)
            if weighted_ec_loss:
                row["ec_loss"] = weighted_ec_loss / max(seen, 1)
            history.append(row)
            if progress_callback is not None:
                progress_callback(row)
        return history

    def state_array(self):
        return self.state.detach().cpu().numpy()

    def write(self, output_prefix, diagnostics=None):
        """Write the shared state, decoder parameters, and identifiers."""
        output_prefix = Path(output_prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(f"{output_prefix}glm_rows.txt", self.data.barcodes, fmt="%s")
        np.savetxt(f"{output_prefix}glm_cols.txt", self.data.transcripts, fmt="%s")
        np.savez_compressed(
            f"{output_prefix}glm_factors.npz",
            left=self.isoform_loadings.detach().cpu().numpy(),
            right=self.state_array(),
        )
        np.savez_compressed(
            f"{output_prefix}hierarchical_parameters.npz",
            gene_intercept=self.gene_intercept.detach().cpu().numpy(),
            gene_loadings=self.gene_loadings.detach().cpu().numpy(),
            isoform_intercept=self.isoform_intercept.detach().cpu().numpy(),
            isoform_loadings=self.isoform_loadings.detach().cpu().numpy(),
            transcript_gene=self.data.transcript_gene,
            genes=self.data.genes,
        )
        payload = {
            "method": "hierarchical",
            "rank": self.rank,
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "n_transcripts": self.n_transcripts,
            **(diagnostics or {}),
        }
        Path(f"{output_prefix}glm_manifest.json").write_text(
            json.dumps({"method": "hierarchical", "diagnostics": payload}, indent=2)
            + "\n"
        )
