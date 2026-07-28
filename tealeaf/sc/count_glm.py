"""Count-aware flat transcript GLM for paired single-cell EC observations."""

from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from tealeaf.sc.hierarchical import _event_edges, _torch


class CountAwareGLM:
    """Flat log-linear transcript decoder with NB gene and paired EC losses."""

    def __init__(
        self,
        data,
        rank,
        selected_genes,
        *,
        gene_family="negative_binomial",
        concentration=10.0,
        target_sum=10_000.0,
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
        self.selected_genes = np.asarray(selected_genes, dtype=np.int64)
        if (
            self.selected_genes.ndim != 1
            or len(self.selected_genes) == 0
            or self.selected_genes.min() < 0
            or self.selected_genes.max() >= self.n_genes
        ):
            raise ValueError("selected_genes are invalid")
        self.selected_genes_t = torch.as_tensor(
            self.selected_genes, dtype=torch.int64, device=self.device
        )
        if gene_family not in {
            "negative_binomial",
            "poisson",
            "standardized_log",
        }:
            raise ValueError(
                "gene_family must be negative_binomial, poisson, "
                "or standardized_log"
            )
        self.gene_family = gene_family
        self.target_sum = float(target_sum)
        if not np.isfinite(self.target_sum) or self.target_sum <= 0:
            raise ValueError("target_sum must be positive and finite")
        concentration = float(concentration)
        if not np.isfinite(concentration) or concentration <= 0:
            raise ValueError("concentration must be positive and finite")
        self.concentration = concentration
        self.gene_library_totals = (
            np.asarray(data.gene_counts.sum(axis=1)).ravel().astype(np.float32)
        )
        if self.gene_family == "standardized_log":
            transformed = data.gene_counts[:, self.selected_genes].astype(np.float32)
            scales = np.zeros_like(self.gene_library_totals)
            active = self.gene_library_totals > 0
            scales[active] = self.target_sum / self.gene_library_totals[active]
            transformed.data *= np.repeat(scales, np.diff(transformed.indptr))
            np.log1p(transformed.data, out=transformed.data)
            means = np.asarray(transformed.mean(axis=0)).ravel()
            second = np.asarray(transformed.power(2).mean(axis=0)).ravel()
            feature_scales = np.sqrt(np.maximum(second - np.square(means), 0.0))
            feature_scales[~np.isfinite(feature_scales) | (feature_scales == 0)] = 1.0
            self.log_gene_means = torch.as_tensor(
                means, dtype=torch.float32, device=self.device
            )
            self.log_gene_scales = torch.as_tensor(
                feature_scales, dtype=torch.float32, device=self.device
            )

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
        scale = initial_state.std(axis=0, keepdims=True)
        scale[scale == 0] = 1
        initial_state /= scale

        gene_totals = np.asarray(data.gene_counts.sum(axis=0)).ravel()
        transcript_intercept = np.log(gene_totals[data.transcript_gene] + 0.5) - np.log(
            np.bincount(data.transcript_gene)[data.transcript_gene]
        )
        transcript_intercept -= transcript_intercept.mean()
        self.state = torch.nn.Parameter(
            torch.as_tensor(initial_state, device=self.device)
        )
        self.transcript_intercept = torch.nn.Parameter(
            torch.as_tensor(
                transcript_intercept, dtype=torch.float32, device=self.device
            )
        )
        self.transcript_loadings = torch.nn.Parameter(
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

    def parameters(self, stage):
        if stage == "gene" or stage == "joint":
            return [
                self.state,
                self.transcript_intercept,
                self.transcript_loadings,
            ]
        if stage == "ec":
            return [self.transcript_intercept, self.transcript_loadings]
        raise ValueError("stage must be gene, ec, or joint")

    def _abundance(self, cell_indices):
        logits = (
            self.transcript_intercept[:, None]
            + self.transcript_loadings @ self.state[cell_indices].T
        ).clamp(-20, 20)
        theta = self.torch.exp(logits.index_select(0, self.transcript_order))
        return theta

    def gene_loss(self, cell_indices):
        torch = self.torch
        cell_indices = np.asarray(cell_indices, dtype=np.int64)
        tensor_cells = torch.as_tensor(
            cell_indices, dtype=torch.int64, device=self.device
        )
        theta = self._abundance(tensor_cells)
        all_gene_mass = torch.segment_reduce(
            theta, reduce="sum", lengths=self.gene_lengths, axis=0
        )
        gene_mass = all_gene_mass.index_select(0, self.selected_genes_t)
        observed = torch.as_tensor(
            self.data.gene_counts[cell_indices][:, self.selected_genes].toarray().T,
            dtype=torch.float32,
            device=self.device,
        )
        totals = observed.sum(dim=0)
        if self.gene_family == "standardized_log":
            library_totals = torch.as_tensor(
                self.gene_library_totals[cell_indices],
                dtype=torch.float32,
                device=self.device,
            )
            observed = torch.log1p(
                observed * (self.target_sum / library_totals.clamp_min(1.0))[None, :]
            )
            predicted = torch.log1p(
                self.target_sum
                * gene_mass
                / all_gene_mass.sum(dim=0, keepdim=True).clamp_min(1e-30)
            )
            observed = (observed - self.log_gene_means[:, None]) / self.log_gene_scales[
                :, None
            ]
            predicted = (
                predicted - self.log_gene_means[:, None]
            ) / self.log_gene_scales[:, None]
            return 0.5 * (observed - predicted).square().mean()
        expected = (
            gene_mass
            / gene_mass.sum(dim=0, keepdim=True).clamp_min(1e-30)
            * totals[None, :]
        ).clamp_min(1e-8)
        if self.gene_family == "poisson":
            negative_log_probability = (
                expected - observed * torch.log(expected) + torch.lgamma(observed + 1)
            )
            return negative_log_probability.sum() / totals.sum().clamp_min(1.0)
        concentration = torch.as_tensor(
            self.concentration, dtype=torch.float32, device=self.device
        )
        log_probability = (
            torch.lgamma(observed + concentration)
            - torch.lgamma(concentration)
            - torch.lgamma(observed + 1)
            + concentration
            * (torch.log(concentration) - torch.log(concentration + expected))
            + observed * (torch.log(expected) - torch.log(concentration + expected))
        )
        return -log_probability.sum() / totals.sum().clamp_min(1.0)

    def ec_loss(self, cell_indices, loading_penalty=0.0):
        torch = self.torch
        cell_indices = np.asarray(cell_indices, dtype=np.int64)
        tensor_cells = torch.as_tensor(
            cell_indices, dtype=torch.int64, device=self.device
        )
        theta = self._abundance(tensor_cells)
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
            if not len(event_counts):
                continue
            edge_events_t = torch.as_tensor(
                edge_events, dtype=torch.int64, device=self.device
            )
            edge_cells_t = torch.as_tensor(
                edge_cells, dtype=torch.int64, device=self.device
            )
            sorted_transcripts = self.original_to_sorted.index_select(
                0,
                torch.as_tensor(
                    edge_transcripts, dtype=torch.int64, device=self.device
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
            values = torch.as_tensor(
                event_counts, dtype=torch.float32, device=self.device
            )
            event_cells_t = torch.as_tensor(
                event_cells, dtype=torch.int64, device=self.device
            )
            total_loss += (
                values
                * (
                    torch.log(normalizer[event_cells_t])
                    - torch.log(predictions.clamp_min(1e-30))
                )
            ).sum()
        loss = total_loss / max(total_molecules, 1.0)
        if loading_penalty:
            loss = loss + 0.5 * float(loading_penalty) * (
                self.transcript_loadings.square().mean()
            )
        return loss

    def center_identifiable_parameters(self):
        with self.torch.no_grad():
            self.transcript_intercept.sub_(self.transcript_intercept.mean())

    def fit_stage(
        self,
        stage,
        *,
        epochs,
        batch_cells,
        learning_rate,
        gene_weight=1.0,
        loading_penalty=0.0,
        seed=0,
        progress_callback=None,
    ):
        torch = self.torch
        optimizer = torch.optim.Adam(self.parameters(stage), lr=float(learning_rate))
        rng = np.random.default_rng(seed)
        history = []
        for epoch in range(int(epochs)):
            started = time.perf_counter()
            order = rng.permutation(self.n_cells)
            totals = {"loss": 0.0, "gene_loss": 0.0, "ec_loss": 0.0}
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
                    ec_component = self.ec_loss(batch, loading_penalty)
                    loss = ec_component
                    if stage == "joint" and gene_weight:
                        gene_component = self.gene_loss(batch)
                        loss = loss + float(gene_weight) * gene_component
                if not bool(torch.isfinite(loss).item()):
                    raise FloatingPointError(
                        f"non-finite {stage} loss at epoch {epoch + 1}"
                    )
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.parameters(stage), max_norm=10.0)
                optimizer.step()
                weight = len(batch)
                totals["loss"] += float(loss.detach().cpu()) * weight
                if gene_component is not None:
                    totals["gene_loss"] += float(gene_component.detach().cpu()) * weight
                if ec_component is not None:
                    totals["ec_loss"] += float(ec_component.detach().cpu()) * weight
                seen += weight
            self.center_identifiable_parameters()
            row = {
                "stage": stage,
                "epoch": epoch + 1,
                "loss": totals["loss"] / max(seen, 1),
                "seconds": time.perf_counter() - started,
            }
            for name in ("gene_loss", "ec_loss"):
                if totals[name]:
                    row[name] = totals[name] / max(seen, 1)
            history.append(row)
            if progress_callback is not None:
                progress_callback(row)
        return history

    def write(self, output_prefix, diagnostics=None):
        output_prefix = Path(output_prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(f"{output_prefix}glm_rows.txt", self.data.barcodes, fmt="%s")
        np.savetxt(f"{output_prefix}glm_cols.txt", self.data.transcripts, fmt="%s")
        np.savez_compressed(
            f"{output_prefix}glm_factors.npz",
            left=self.transcript_loadings.detach().cpu().numpy(),
            right=self.state.detach().cpu().numpy(),
        )
        np.savez_compressed(
            f"{output_prefix}count_glm_parameters.npz",
            transcript_intercept=self.transcript_intercept.detach().cpu().numpy(),
            transcript_loadings=self.transcript_loadings.detach().cpu().numpy(),
            selected_genes=self.selected_genes,
            genes=self.data.genes,
            transcript_gene=self.data.transcript_gene,
        )
        payload = {
            "method": "count_aware_flat_glm",
            "rank": self.rank,
            "gene_family": self.gene_family,
            "concentration": self.concentration,
            "target_sum": self.target_sum,
            "n_cells": self.n_cells,
            "n_selected_genes": int(len(self.selected_genes)),
            **(diagnostics or {}),
        }
        Path(f"{output_prefix}glm_manifest.json").write_text(
            json.dumps({"method": payload["method"], "diagnostics": payload}, indent=2)
            + "\n"
        )
