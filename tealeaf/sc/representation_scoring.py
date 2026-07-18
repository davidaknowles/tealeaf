"""Benchmarks for fitted cell representations against reference labels."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_glm_factors(prefix):
    """Load cell IDs, transcript IDs, and both factors from a GLM fit."""
    prefix = Path(prefix)
    cell_ids = np.loadtxt(f"{prefix}glm_rows.txt", dtype=str, ndmin=1)
    transcript_ids = np.loadtxt(f"{prefix}glm_cols.txt", dtype=str, ndmin=1)
    with np.load(f"{prefix}glm_factors.npz") as factors:
        left = factors["left"]
        right = factors["right"]
    if right.shape[0] != len(cell_ids):
        raise ValueError(
            f"factor rows ({right.shape[0]}) do not match cell IDs ({len(cell_ids)})"
        )
    if left.shape[0] != len(transcript_ids):
        raise ValueError(
            f"factor columns ({left.shape[0]}) do not match transcripts "
            f"({len(transcript_ids)})"
        )
    if left.shape[1] != right.shape[1]:
        raise ValueError("left and right factors have different ranks")
    return cell_ids, transcript_ids, left, right


def load_factor_representation(prefix):
    """Load cell identities and the right factor from a single-cell GLM fit."""
    prefix = Path(prefix)
    cell_ids = np.loadtxt(f"{prefix}glm_rows.txt", dtype=str, ndmin=1)
    with np.load(f"{prefix}glm_factors.npz") as factors:
        right = factors["right"]
    if right.shape[0] != len(cell_ids):
        raise ValueError(
            f"factor rows ({right.shape[0]}) do not match cell IDs ({len(cell_ids)})"
        )
    return cell_ids, right


def align_reference_metadata(cell_ids, labels, groups=None):
    """Return fitted row positions aligned to unique reference metadata."""
    labels = labels.loc[:, ["cell_id", "label"]].dropna().drop_duplicates()
    conflicting = labels.groupby("cell_id")["label"].nunique()
    if (conflicting > 1).any():
        raise ValueError("reference labels contain conflicting labels for a cell ID")
    labels = labels.drop_duplicates("cell_id").set_index("cell_id")
    if pd.Index(cell_ids).has_duplicates:
        raise ValueError("fitted cell IDs must be unique")
    positions = pd.Series(
        np.arange(len(cell_ids)), index=pd.Index(cell_ids, name="cell_id")
    )
    aligned = labels.join(positions.rename("position"), how="inner")
    if groups is not None:
        groups = groups.loc[:, ["cell_id", "group"]].dropna().drop_duplicates()
        group_conflicts = groups.groupby("cell_id")["group"].nunique()
        if (group_conflicts > 1).any():
            raise ValueError("reference groups contain conflicting groups for a cell ID")
        aligned = aligned.join(
            groups.drop_duplicates("cell_id").set_index("cell_id"), how="inner"
        )
    return (
        aligned["position"].to_numpy(dtype=int),
        aligned["label"].to_numpy(dtype=str),
        None if groups is None else aligned["group"].to_numpy(dtype=str),
        aligned.index.to_numpy(dtype=str),
    )


def align_reference_labels(cell_ids, representation, labels, groups=None):
    """Align a representation to unique reference labels by cell ID."""
    positions, aligned_labels, aligned_groups, aligned_ids = (
        align_reference_metadata(cell_ids, labels, groups)
    )
    return (
        representation[positions],
        aligned_labels,
        aligned_groups,
        aligned_ids,
    )


def aggregate_transcript_loadings(
    left, transcript_ids, transcript_to_gene, eligible_genes
):
    """Sum transcript loadings into an ordered eligible-gene matrix."""
    mapping = transcript_to_gene.loc[:, ["transcript_id", "gene_id"]].dropna()
    conflicts = mapping.groupby("transcript_id")["gene_id"].nunique()
    if (conflicts > 1).any():
        raise ValueError("transcript-to-gene mapping assigns transcripts to multiple genes")
    transcript_lookup = mapping.drop_duplicates("transcript_id").set_index(
        "transcript_id"
    )["gene_id"]
    mapped_genes = transcript_lookup.reindex(np.asarray(transcript_ids, dtype=str))
    eligible_genes = np.asarray(eligible_genes, dtype=str)
    if pd.Index(eligible_genes).has_duplicates:
        raise ValueError("eligible gene IDs must be unique")
    gene_lookup = pd.Series(np.arange(len(eligible_genes)), index=eligible_genes)
    gene_indices = gene_lookup.reindex(mapped_genes.to_numpy()).to_numpy()
    keep = ~pd.isna(gene_indices)
    gene_indices = gene_indices[keep].astype(int)
    gene_left = np.zeros(
        (len(eligible_genes), left.shape[1]), dtype=np.asarray(left).dtype
    )
    np.add.at(gene_left, gene_indices, np.asarray(left)[keep])
    diagnostics = {
        "n_transcripts": int(len(transcript_ids)),
        "n_mapped_transcripts": int(keep.sum()),
        "n_eligible_genes": int(len(eligible_genes)),
        "n_genes_with_loadings": int(np.count_nonzero(np.linalg.norm(gene_left, axis=1))),
    }
    return gene_left, diagnostics


class FactorizedGeneExpression:
    """Stream normalized log-gene abundance from compact cell/gene factors."""

    def __init__(self, right, gene_left, *, target_sum=10_000.0,
                 batch_cells=2048, device="auto"):
        import torch

        self.torch = torch
        self.right = np.asarray(right, dtype=np.float32)
        self.gene_left = np.asarray(gene_left, dtype=np.float32)
        if self.right.ndim != 2 or self.gene_left.ndim != 2:
            raise ValueError("right and gene_left must be matrices")
        if self.right.shape[1] != self.gene_left.shape[1]:
            raise ValueError("right and gene_left have different ranks")
        self.target_sum = float(target_sum)
        self.batch_cells = max(1, int(batch_cells))
        requested = "cuda" if device == "auto" and torch.cuda.is_available() else device
        requested = "cpu" if requested == "auto" else requested
        self.device = torch.device(requested)
        if self.device.type == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA reconstruction was requested but is unavailable")
        self._gene_left = torch.as_tensor(self.gene_left, device=self.device)
        total_loading = self.gene_left.sum(axis=0, dtype=np.float64)
        self.row_totals = self.right.astype(np.float64) @ total_loading
        self.active = np.isfinite(self.row_totals) & (self.row_totals > 0)

    def _bounds(self, n_rows, minimum_last=1):
        starts = list(range(0, n_rows, self.batch_cells))
        if len(starts) > 1 and n_rows - starts[-1] < minimum_last:
            starts.pop()
        return [(start, starts[i + 1] if i + 1 < len(starts) else n_rows)
                for i, start in enumerate(starts)]

    def iter_blocks(self, gene_indices=None, minimum_last=1):
        """Yield active-row positions and normalized log1p gene blocks."""
        active_positions = np.flatnonzero(self.active)
        if gene_indices is None:
            selected_left = self._gene_left
        else:
            indices = self.torch.as_tensor(
                np.asarray(gene_indices, dtype=np.int64),
                dtype=self.torch.int64,
                device=self.device,
            )
            selected_left = self._gene_left.index_select(0, indices)
        for start, stop in self._bounds(len(active_positions), minimum_last):
            positions = active_positions[start:stop]
            right_block = self.torch.as_tensor(
                self.right[positions], device=self.device
            )
            abundance = self.torch.relu(right_block @ selected_left.T)
            scales = self.torch.as_tensor(
                self.target_sum / self.row_totals[positions],
                dtype=abundance.dtype,
                device=self.device,
            )
            abundance.mul_(scales[:, None]).log1p_()
            yield positions, abundance.detach().cpu().numpy()

    def mean_variance(self):
        """Compute per-gene log-abundance moments without materializing cells by genes."""
        means = np.zeros(self.gene_left.shape[0], dtype=np.float64)
        centered_ss = np.zeros_like(means)
        n_seen = 0
        n_rows = int(self.active.sum())
        if n_rows < 2:
            raise ValueError("at least two active cells are required")
        for _, block in self.iter_blocks():
            block = np.asarray(block, dtype=np.float64)
            block_n = len(block)
            block_mean = block.mean(axis=0)
            block_ss = np.square(block - block_mean).sum(axis=0)
            delta = block_mean - means
            combined_n = n_seen + block_n
            centered_ss += block_ss + np.square(delta) * n_seen * block_n / combined_n
            means += delta * block_n / combined_n
            n_seen = combined_n
        variances = centered_ss / (n_rows - 1)
        return means, variances

    def fit_pca(self, gene_indices, n_components=30):
        """Fit and transform incremental PCA in repeated reconstruction passes."""
        from sklearn.decomposition import IncrementalPCA

        n_active = int(self.active.sum())
        n_components = min(int(n_components), len(gene_indices), n_active - 1)
        if n_components < 1:
            raise ValueError("PCA requires at least one component")
        pca = IncrementalPCA(
            n_components=n_components, batch_size=self.batch_cells
        )
        for _, block in self.iter_blocks(
            gene_indices, minimum_last=n_components
        ):
            pca.partial_fit(block)
        embedding = np.full(
            (len(self.right), n_components), np.nan, dtype=np.float32
        )
        for positions, block in self.iter_blocks(gene_indices):
            embedding[positions] = pca.transform(block).astype(np.float32)
        return embedding, pca


def log_gene_pca_embedding(
    right,
    gene_left,
    gene_ids,
    *,
    target_sum=10_000.0,
    n_hvg=2_000,
    n_components=30,
    batch_cells=2048,
    device="auto",
):
    """Build a per-fit HVG PCA embedding from normalized log gene abundance."""
    expression = FactorizedGeneExpression(
        right,
        gene_left,
        target_sum=target_sum,
        batch_cells=batch_cells,
        device=device,
    )
    means, variances = expression.mean_variance()
    total_variance = float(variances.sum(dtype=np.float64))
    mean_square = float(np.square(means).sum(dtype=np.float64))
    relative_variance = total_variance / max(mean_square, np.finfo(float).tiny)
    if not np.isfinite(relative_variance) or relative_variance <= 1e-8:
        raise ValueError(
            "normalized log-gene representation is collapsed: relative "
            f"between-cell variance={relative_variance:.3g}"
        )
    informative = np.isfinite(variances) & (variances > 0)
    candidates = np.flatnonzero(informative)
    if candidates.size == 0:
        raise ValueError("no variable fitted genes are available for PCA")
    n_hvg = min(int(n_hvg), len(candidates))
    selected = candidates[np.argsort(variances[candidates])[-n_hvg:][::-1]]
    embedding, pca = expression.fit_pca(selected, n_components=n_components)
    selected_variance = float(variances[selected].sum(dtype=np.float64))
    explained_variance = float(
        np.asarray(pca.explained_variance_, dtype=np.float64).sum()
        / selected_variance
    )
    diagnostics = {
        "representation": "log1p_gene_pca",
        "target_sum": float(target_sum),
        "n_input_cells": int(len(right)),
        "n_active_cells": int(expression.active.sum()),
        "n_hvg": int(n_hvg),
        "pca_components": int(embedding.shape[1]),
        "pca_explained_variance": float(np.clip(explained_variance, 0.0, 1.0)),
        "total_gene_variance": total_variance,
        "relative_gene_variance": relative_variance,
        "reconstruction_device": str(expression.device),
        "mean_library_total": float(np.mean(expression.row_totals[expression.active])),
    }
    return embedding, np.asarray(gene_ids, dtype=str)[selected], expression.active, diagnostics


def score_embedding(
    embedding,
    labels,
    groups=None,
    *,
    valid_mask=None,
    n_splits=5,
    random_state=0,
    silhouette_sample_size=10_000,
):
    """Score grouped label transfer and clustering in a supplied embedding."""
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        accuracy_score,
        adjusted_rand_score,
        balanced_accuracy_score,
        f1_score,
        normalized_mutual_info_score,
        silhouette_score,
    )
    from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    embedding = np.asarray(embedding, dtype=np.float32)
    labels = np.asarray(labels, dtype=str)
    keep = np.isfinite(embedding).all(axis=1)
    if valid_mask is not None:
        keep &= np.asarray(valid_mask, dtype=bool)
    coverage = float(keep.mean()) if len(keep) else 0.0
    embedding = embedding[keep]
    labels = labels[keep]
    if groups is not None:
        groups = np.asarray(groups, dtype=str)[keep]
    if embedding.shape[0] == 0:
        raise ValueError("embedding has no finite, valid labeled cells")

    class_counts = pd.Series(labels).value_counts()
    splits = min(int(n_splits), int(class_counts.min()))
    if groups is not None:
        splits = min(splits, len(np.unique(groups)))
        splitter = StratifiedGroupKFold(
            n_splits=splits, shuffle=True, random_state=random_state
        )
        split_iterator = splitter.split(embedding, labels, groups)
    else:
        splitter = StratifiedKFold(
            n_splits=splits, shuffle=True, random_state=random_state
        )
        split_iterator = splitter.split(embedding, labels)
    if splits < 2:
        raise ValueError("at least two folds are required for label scoring")

    fold_rows = []
    for fold, (train, test) in enumerate(split_iterator):
        classifier = make_pipeline(
            StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=1000),
        )
        classifier.fit(embedding[train], labels[train])
        predicted = classifier.predict(embedding[test])
        fold_rows.append(
            {
                "fold": fold,
                "n_train": len(train),
                "n_test": len(test),
                "accuracy": accuracy_score(labels[test], predicted),
                "balanced_accuracy": balanced_accuracy_score(labels[test], predicted),
                "macro_f1": f1_score(labels[test], predicted, average="macro"),
            }
        )
    folds = pd.DataFrame(fold_rows)

    clusters = MiniBatchKMeans(
        n_clusters=len(np.unique(labels)),
        batch_size=4096,
        n_init=10,
        random_state=random_state,
    ).fit_predict(embedding)
    silhouette_size = min(int(silhouette_sample_size), len(labels))
    report = {
        "n_reference_cells": int(len(keep)),
        "n_scored_cells": int(keep.sum()),
        "active_finite_coverage": coverage,
        "n_labels": int(len(np.unique(labels))),
        "n_groups": None if groups is None else int(len(np.unique(groups))),
        "cv_folds": int(splits),
        "embedding_dimensions": int(embedding.shape[1]),
        "silhouette_sample_size": int(silhouette_size),
        "accuracy_mean": float(folds["accuracy"].mean()),
        "accuracy_sd": float(folds["accuracy"].std(ddof=1)),
        "balanced_accuracy_mean": float(folds["balanced_accuracy"].mean()),
        "balanced_accuracy_sd": float(folds["balanced_accuracy"].std(ddof=1)),
        "macro_f1_mean": float(folds["macro_f1"].mean()),
        "macro_f1_sd": float(folds["macro_f1"].std(ddof=1)),
        "adjusted_rand_index": float(adjusted_rand_score(labels, clusters)),
        "normalized_mutual_info": float(normalized_mutual_info_score(labels, clusters)),
        "reference_label_silhouette": float(
            silhouette_score(
                embedding,
                labels,
                sample_size=silhouette_size,
                random_state=random_state,
            )
        ),
        "kmeans_cluster_silhouette": float(
            silhouette_score(
                embedding,
                clusters,
                sample_size=silhouette_size,
                random_state=random_state,
            )
        ),
    }
    return report, folds


def score_representation(
    representation,
    labels,
    groups=None,
    *,
    n_splits=5,
    random_state=0,
    pca_components=30,
    silhouette_sample_size=10_000,
):
    """PCA-transform and score a generic representation."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    representation = np.asarray(representation, dtype=np.float32)
    finite = np.isfinite(representation).all(axis=1)
    active = np.linalg.norm(np.nan_to_num(representation), axis=1) > 0
    valid = finite & active
    n_components = min(
        int(pca_components), representation.shape[1], int(valid.sum()) - 1
    )
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=random_state)
    embedding = np.full((len(representation), n_components), np.nan, dtype=np.float32)
    embedding[valid] = pca.fit_transform(
        StandardScaler().fit_transform(representation[valid])
    )
    report, folds = score_embedding(
        embedding,
        labels,
        groups,
        valid_mask=valid,
        n_splits=n_splits,
        random_state=random_state,
        silhouette_sample_size=silhouette_sample_size,
    )
    report.update(
        pca_components=int(n_components),
        pca_explained_variance=float(pca.explained_variance_ratio_.sum()),
    )
    return report, folds


def write_embedding(embedding, cell_ids, genes, diagnostics, output_prefix):
    """Write PCA coordinates and preprocessing metadata."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(f"{output_prefix}log_gene_pca.npz", embedding=embedding)
    np.savetxt(f"{output_prefix}log_gene_pca_cells.txt", cell_ids, fmt="%s")
    np.savetxt(f"{output_prefix}log_gene_pca_genes.txt", genes, fmt="%s")
    Path(f"{output_prefix}log_gene_pca.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n"
    )


def write_score(report, folds, output_prefix):
    """Write a JSON summary and per-fold CSV."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{output_prefix}label_scores.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    folds.to_csv(f"{output_prefix}label_score_folds.csv", index=False)
