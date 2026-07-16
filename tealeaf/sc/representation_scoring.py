"""Benchmarks for cell representations against reference cell labels."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_factor_representation(prefix):
    """Load cell identities and the right factor from a single-cell GLM fit."""
    prefix = Path(prefix)
    rows = np.loadtxt(f"{prefix}glm_rows.txt", dtype=str, ndmin=1)
    with np.load(f"{prefix}glm_factors.npz") as factors:
        representation = factors["right"]
    if representation.shape[0] != len(rows):
        raise ValueError(
            f"factor rows ({representation.shape[0]}) do not match cell IDs ({len(rows)})"
        )
    return rows, representation


def align_reference_labels(cell_ids, representation, labels, groups=None):
    """Align a representation to unique reference labels by cell ID."""
    labels = labels.loc[:, ["cell_id", "label"]].dropna().drop_duplicates()
    conflicting = labels.groupby("cell_id")["label"].nunique()
    if (conflicting > 1).any():
        raise ValueError("reference labels contain conflicting labels for a cell ID")
    labels = labels.drop_duplicates("cell_id").set_index("cell_id")
    positions = pd.Series(np.arange(len(cell_ids)), index=pd.Index(cell_ids, name="cell_id"))
    aligned = labels.join(positions.rename("position"), how="inner")
    if groups is not None:
        groups = groups.loc[:, ["cell_id", "group"]].dropna().drop_duplicates()
        group_conflicts = groups.groupby("cell_id")["group"].nunique()
        if (group_conflicts > 1).any():
            raise ValueError("reference groups contain conflicting groups for a cell ID")
        aligned = aligned.join(
            groups.drop_duplicates("cell_id").set_index("cell_id"), how="inner"
        )
    positions_array = aligned["position"].to_numpy(dtype=int)
    return (
        representation[positions_array],
        aligned["label"].to_numpy(dtype=str),
        None if groups is None else aligned["group"].to_numpy(dtype=str),
        aligned.index.to_numpy(dtype=str),
    )


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
    """Score supervised label transfer and unsupervised label recovery.

    Supervised metrics use stratified held-out folds. When groups are supplied,
    no group is shared between training and test data. Unsupervised metrics use
    mini-batch k-means with the reference number of classes.
    """
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
    from sklearn.decomposition import PCA
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    representation = np.asarray(representation, dtype=np.float32)
    labels = np.asarray(labels, dtype=str)
    finite = np.isfinite(representation).all(axis=1)
    active = np.linalg.norm(np.nan_to_num(representation), axis=1) > 0
    keep = finite & active
    coverage = float(keep.mean()) if len(keep) else 0.0
    representation = representation[keep]
    labels = labels[keep]
    if groups is not None:
        groups = np.asarray(groups, dtype=str)[keep]
    if representation.shape[0] == 0:
        raise ValueError("representation has no finite, nonzero labeled cells")

    class_counts = pd.Series(labels).value_counts()
    splits = min(int(n_splits), int(class_counts.min()))
    if groups is not None:
        splits = min(splits, len(np.unique(groups)))
        splitter = StratifiedGroupKFold(
            n_splits=splits, shuffle=True, random_state=random_state
        )
        split_iterator = splitter.split(representation, labels, groups)
    else:
        splitter = StratifiedKFold(
            n_splits=splits, shuffle=True, random_state=random_state
        )
        split_iterator = splitter.split(representation, labels)
    if splits < 2:
        raise ValueError("at least two folds are required for label scoring")

    fold_rows = []
    for fold, (train, test) in enumerate(split_iterator):
        classifier = make_pipeline(
            StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=1000),
        )
        classifier.fit(representation[train], labels[train])
        predicted = classifier.predict(representation[test])
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

    scaled_representation = StandardScaler().fit_transform(representation)
    n_components = min(
        int(pca_components), representation.shape[1], representation.shape[0] - 1
    )
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=random_state)
    pca_representation = pca.fit_transform(scaled_representation)
    clusters = MiniBatchKMeans(
        n_clusters=len(np.unique(labels)),
        batch_size=4096,
        n_init=10,
        random_state=random_state,
    ).fit_predict(pca_representation)
    silhouette_size = min(int(silhouette_sample_size), len(labels))
    label_silhouette = silhouette_score(
        pca_representation,
        labels,
        sample_size=silhouette_size,
        random_state=random_state,
    )
    cluster_silhouette = silhouette_score(
        pca_representation,
        clusters,
        sample_size=silhouette_size,
        random_state=random_state,
    )
    report = {
        "n_reference_cells": int(len(keep)),
        "n_scored_cells": int(keep.sum()),
        "active_finite_coverage": coverage,
        "n_labels": int(len(np.unique(labels))),
        "n_groups": None if groups is None else int(len(np.unique(groups))),
        "cv_folds": int(splits),
        "pca_components": int(n_components),
        "pca_explained_variance": float(pca.explained_variance_ratio_.sum()),
        "silhouette_sample_size": int(silhouette_size),
        "accuracy_mean": float(folds["accuracy"].mean()),
        "accuracy_sd": float(folds["accuracy"].std(ddof=1)),
        "balanced_accuracy_mean": float(folds["balanced_accuracy"].mean()),
        "balanced_accuracy_sd": float(folds["balanced_accuracy"].std(ddof=1)),
        "macro_f1_mean": float(folds["macro_f1"].mean()),
        "macro_f1_sd": float(folds["macro_f1"].std(ddof=1)),
        "adjusted_rand_index": float(adjusted_rand_score(labels, clusters)),
        "normalized_mutual_info": float(normalized_mutual_info_score(labels, clusters)),
        "reference_label_silhouette": float(label_silhouette),
        "kmeans_cluster_silhouette": float(cluster_silhouette),
    }
    return report, folds


def write_score(report, folds, output_prefix):
    """Write a JSON summary and per-fold CSV."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{output_prefix}label_scores.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    folds.to_csv(f"{output_prefix}label_score_folds.csv", index=False)
