"""Reusable data preparation and count-fold CV for scalable single-cell GLMs."""

from __future__ import annotations

from dataclasses import dataclass
import gc
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

from tealeaf.sc import glm_solvers, sc_utils


@dataclass
class PreparedGLMData:
    counts: sp.csr_matrix
    compatibility: sp.csr_matrix
    barcodes: np.ndarray
    features: np.ndarray


def prepare_alevin_glm_data(
    alevin_dir,
    salmon_ref,
    *,
    ec_design="legacy",
    regularization_target="phi",
    min_eq=5,
    probability_file=None,
    weight_cache=None,
):
    """Load and filter the matrices shared by single-cell fitting and CV."""
    alevin_dir = Path(alevin_dir)
    map_cache = alevin_dir / "gene_eqclass.npz"
    if map_cache.is_file():
        ec_transcript = sp.load_npz(map_cache).tocsr()
    else:
        _, _, ecs = sc_utils.read_alevin_ec(alevin_dir / "gene_eqclass.txt.gz")
        ec_transcript = sc_utils.to_coo(
            [ecs[index] for index in range(len(ecs))]
        ).tocsr()
        sp.save_npz(map_cache, ec_transcript)

    count_cache = alevin_dir / "geqc_counts.npz"
    if count_cache.is_file():
        counts = sp.load_npz(count_cache).tocsr()
    else:
        counts = scipy.io.mmread(alevin_dir / "geqc_counts.mtx").tocsr()
        sp.save_npz(count_cache, counts)
    with open(alevin_dir / "quants_mat_cols.txt") as handle:
        features = np.asarray([line.strip() for line in handle])
    with open(alevin_dir / "quants_mat_rows.txt") as handle:
        barcodes = np.asarray([line.strip() for line in handle])

    transcript_lengths = sc_utils.get_transcript_lengths(Path(salmon_ref))
    _, full_weights = sc_utils.get_feature_weights(features, transcript_lengths)
    probability_file = probability_file or (
        alevin_dir / "gene_eqclass_probs.tsv.gz"
    )
    weight_cache = weight_cache or (
        alevin_dir / "gene_eqclass_fixed_weights.npz"
    )
    phi_design = sc_utils.glm_design_matrix(
        ec_transcript,
        full_weights,
        parameterization="phi",
        design=ec_design,
        probability_file=probability_file,
        cache_file=weight_cache,
    )

    aggregate = sc_utils.sparse_sum(counts, 0)
    transcript_count = aggregate @ ec_transcript
    ec_keep = aggregate >= float(min_eq)
    feature_keep = transcript_count > 0
    filtered_counts = counts[:, ec_keep].tocsr()
    filtered_phi_design = phi_design[ec_keep, :][:, feature_keep]
    filtered_features = features[feature_keep]
    _, weights = sc_utils.get_feature_weights(
        filtered_features, transcript_lengths
    )
    compatibility = sc_utils.parameterize_glm_design(
        filtered_phi_design,
        weights,
        regularization_target,
        normalize_columns=ec_design in {"binary", "weighted"},
    ).tocsr()
    return PreparedGLMData(
        filtered_counts, compatibility, barcodes, filtered_features
    )


def sample_nonempty_cells(counts, n_cells, seed=0):
    """Select a reproducible random subset of nonempty cell rows."""
    totals = np.asarray(counts.sum(axis=1)).ravel()
    eligible = np.flatnonzero(totals > 0)
    if n_cells is None or int(n_cells) >= len(eligible):
        return eligible
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(eligible, size=int(n_cells), replace=False))


def split_count_folds(counts, n_folds=3, seed=0):
    """Randomly partition integer molecule counts into sparse validation folds."""
    if int(n_folds) < 2:
        raise ValueError("n_folds must be at least two")
    counts = counts.tocsr()
    integer_counts = np.rint(counts.data).astype(np.int64)
    if np.any(integer_counts < 0) or not np.allclose(counts.data, integer_counts):
        raise ValueError("count-fold CV requires nonnegative integer counts")
    rng = np.random.default_rng(seed)
    remaining = integer_counts.copy()
    fold_data = []
    for fold in range(int(n_folds) - 1):
        draw = rng.binomial(remaining, 1.0 / (int(n_folds) - fold))
        fold_data.append(draw)
        remaining -= draw
    fold_data.append(remaining)
    folds = []
    for data in fold_data:
        fold = sp.csr_matrix(
            (data, counts.indices.copy(), counts.indptr.copy()),
            shape=counts.shape,
        )
        fold.eliminate_zeros()
        folds.append(fold)
    return folds


def _b_transpose_times(data, vectors):
    result = data.torch.zeros(
        (data.n_cells, vectors.shape[1]), device=data.device
    )
    a_vectors = data.a_times(vectors)
    for start, stop, block, _ in data.blocks():
        result[start:stop] = data.torch.sparse.mm(block, a_vectors)
    return result


def estimate_b_spectral(data, power_iter=10, seed=0):
    """Estimate top singular vectors of B=A.T@C without materializing B."""
    torch = data.torch
    generator = torch.Generator(device=data.device)
    generator.manual_seed(seed)
    right = torch.randn(
        (data.n_cells, 1), generator=generator, device=data.device
    )
    right /= torch.linalg.vector_norm(right).clamp_min(1e-12)
    for _ in range(int(power_iter)):
        left = data.b_times(right)
        left /= torch.linalg.vector_norm(left).clamp_min(1e-12)
        right = _b_transpose_times(data, left)
        right /= torch.linalg.vector_norm(right).clamp_min(1e-12)
    image = data.b_times(right)
    singular = torch.linalg.vector_norm(image).clamp_min(1e-12)
    left = image / singular
    return float(singular.item()), left, right


def hyperparameter_scale(counts, compatibility, method, *, device="auto",
                         batch_cells=4096, power_iter=10, seed=0):
    """Return lambda_max or a rank-one line-search tau reference."""
    data = glm_solvers.SparseGLM(
        counts, compatibility, device=device, batch_cells=batch_cells
    )
    singular, left, _ = estimate_b_spectral(data, power_iter, seed)
    if method == "admm_factorized":
        return singular
    if method == "frank_wolfe_penalized":
        prediction = data.a_times(left)
        curvature = float((prediction * prediction).sum().item())
        return singular / max(curvature, 1e-12)
    raise ValueError(f"CV scale is not defined for method {method}")


def _open_boundary_direction(method, multipliers, best_multiplier):
    ordered = sorted(float(value) for value in multipliers)
    if len(ordered) < 2:
        return None
    lower_needs_expansion = best_multiplier == ordered[0] and ordered[0] != 0.0
    upper_is_lambda_max = method == "admm_factorized" and ordered[-1] >= 1.0
    upper_needs_expansion = (
        best_multiplier == ordered[-1] and not upper_is_lambda_max
    )
    if lower_needs_expansion:
        return "lower"
    if upper_needs_expansion:
        return "upper"
    return None


def _best_on_open_boundary(method, multipliers, best_multiplier):
    return _open_boundary_direction(method, multipliers, best_multiplier) is not None


def _expanded_candidate(method, multipliers, direction, factor):
    if float(factor) <= 1:
        raise ValueError("grid expansion factor must be greater than one")
    ordered = sorted(float(value) for value in multipliers)
    if direction == "lower":
        return ordered[0] / float(factor)
    if direction == "upper":
        candidate = ordered[-1] * float(factor)
        if method == "admm_factorized":
            candidate = min(candidate, 1.0)
        return candidate
    raise ValueError(f"invalid grid expansion direction: {direction}")


def cross_validate_glm(
    counts,
    compatibility,
    method,
    multipliers,
    *,
    n_folds=3,
    seed=0,
    device="auto",
    batch_cells=4096,
    power_iter=10,
    fit_kwargs=None,
):
    """Tune a scale-free lambda fraction or tau multiplier by count folds."""
    fit_kwargs = dict(fit_kwargs or {})
    folds = split_count_folds(counts, n_folds=n_folds, seed=seed)
    rows = []
    scales = []
    for fold_index, validation_counts in enumerate(folds):
        training_counts = sum(
            (fold for index, fold in enumerate(folds) if index != fold_index),
            start=sp.csr_matrix(counts.shape),
        ).tocsr()
        scale = hyperparameter_scale(
            training_counts,
            compatibility,
            method,
            device=device,
            batch_cells=batch_cells,
            power_iter=power_iter,
            seed=seed + fold_index,
        )
        scales.append(scale)
        validation = glm_solvers.SparseGLM(
            validation_counts,
            compatibility,
            device=device,
            batch_cells=batch_cells,
        )
        for multiplier in multipliers:
            kwargs = dict(fit_kwargs)
            if method == "admm_factorized":
                kwargs["regularization"] = float(multiplier) * scale
            elif method == "frank_wolfe_penalized":
                kwargs["tau"] = float(multiplier) * scale
            result = glm_solvers.fit_glm(
                training_counts,
                compatibility,
                method,
                device=device,
                batch_cells=batch_cells,
                **kwargs,
            )
            validation_loss = validation.loss_for_factors(
                result.left, result.right
            ) / validation.n_nonempty
            rows.append(
                {
                    "fold": fold_index,
                    "multiplier": float(multiplier),
                    "scale": scale,
                    "value": float(multiplier) * scale,
                    "validation_loss_per_cell": validation_loss,
                    "iterations": result.diagnostics.get("iterations"),
                    "converged": result.diagnostics.get("converged"),
                }
            )
            del result
            gc.collect()
            if device == "cuda" or (
                device == "auto" and glm_solvers._torch().cuda.is_available()
            ):
                glm_solvers._torch().cuda.empty_cache()
    means = {}
    for multiplier in multipliers:
        values = [
            row["validation_loss_per_cell"]
            for row in rows
            if row["multiplier"] == float(multiplier)
        ]
        means[float(multiplier)] = float(np.mean(values))
    best_multiplier = min(means, key=means.get)
    best_on_boundary = _best_on_open_boundary(
        method, multipliers, best_multiplier
    )
    return {
        "method": method,
        "n_folds": int(n_folds),
        "multipliers": [float(value) for value in multipliers],
        "fold_scales": scales,
        "mean_validation_loss": means,
        "best_multiplier": best_multiplier,
        "best_on_boundary": best_on_boundary,
        "fold_results": rows,
    }


def cross_validate_glm_adaptive_grid(
    counts,
    compatibility,
    method,
    multipliers,
    *,
    max_grid_expansions=4,
    grid_expansion_factor=4.0,
    **kwargs,
):
    """Cross-validate and extend an open-boundary grid one candidate at a time."""
    if int(max_grid_expansions) < 0:
        raise ValueError("max_grid_expansions must be nonnegative")
    if float(grid_expansion_factor) <= 1:
        raise ValueError("grid_expansion_factor must be greater than one")
    candidates = sorted(set(float(value) for value in multipliers))
    if len(candidates) < 2 or any(value < 0 for value in candidates):
        raise ValueError("adaptive CV requires at least two nonnegative candidates")

    report = cross_validate_glm(
        counts, compatibility, method, candidates, **kwargs
    )
    history = [
        {
            "round": 0,
            "added_multiplier": None,
            "best_multiplier": report["best_multiplier"],
            "boundary_direction": _open_boundary_direction(
                method, candidates, report["best_multiplier"]
            ),
        }
    ]
    expansions = 0
    while expansions < int(max_grid_expansions):
        direction = _open_boundary_direction(
            method, candidates, report["best_multiplier"]
        )
        if direction is None:
            break
        candidate = _expanded_candidate(
            method, candidates, direction, grid_expansion_factor
        )
        if candidate in candidates:
            break
        extension = cross_validate_glm(
            counts, compatibility, method, [candidate], **kwargs
        )
        candidates.append(candidate)
        candidates.sort()
        report["fold_results"].extend(extension["fold_results"])
        report["mean_validation_loss"][candidate] = extension[
            "mean_validation_loss"
        ][candidate]
        report["best_multiplier"] = min(
            report["mean_validation_loss"],
            key=report["mean_validation_loss"].get,
        )
        report["multipliers"] = candidates.copy()
        report["best_on_boundary"] = _best_on_open_boundary(
            method, candidates, report["best_multiplier"]
        )
        expansions += 1
        history.append(
            {
                "round": expansions,
                "added_multiplier": candidate,
                "best_multiplier": report["best_multiplier"],
                "boundary_direction": _open_boundary_direction(
                    method, candidates, report["best_multiplier"]
                ),
            }
        )

    report.update(
        grid_expansions=expansions,
        max_grid_expansions=int(max_grid_expansions),
        grid_expansion_factor=float(grid_expansion_factor),
        grid_exhausted=bool(report["best_on_boundary"]),
        grid_history=history,
    )
    return report
