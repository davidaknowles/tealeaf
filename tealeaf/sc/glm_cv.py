"""Reusable data preparation and count-fold CV for scalable single-cell GLMs."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import gc
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

from tealeaf.data.alevin import load_alevin_counts, load_alevin_structure
from tealeaf.sc import glm_solvers, sc_utils


@dataclass
class PreparedGLMData:
    counts: sp.csr_matrix
    compatibility: sp.csr_matrix
    barcodes: np.ndarray
    features: np.ndarray
    cell_umi_totals: np.ndarray
    metadata: dict | None = None
    cv_raw_counts: sp.csr_matrix | None = None


def _row_normalize(matrix):
    matrix = matrix.tocsr().astype(np.float32)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    inverse = np.zeros_like(totals, dtype=np.float32)
    inverse[totals > 0] = 1.0 / totals[totals > 0]
    return (sp.diags(inverse) @ matrix).tocsr(), totals


def _read_primer_pairs(pair_file):
    required = ("cell_id", "polydt_barcode", "ranhex_barcode")
    with open(pair_file, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or any(name not in reader.fieldnames for name in required):
            raise ValueError(
                "primer pair table must contain cell_id, polydt_barcode, and "
                "ranhex_barcode columns"
            )
        rows = []
        pair_to_cell = {}
        half_to_pair = {}
        for row in reader:
            cell_id, polydt, ranhex = (row[name] for name in required)
            pair = (polydt, ranhex)
            if pair in pair_to_cell:
                continue
            for barcode in pair:
                previous = half_to_pair.get(barcode)
                if previous is not None and previous != pair:
                    raise ValueError(
                        "a half-cell barcode occurs in conflicting primer pairs"
                    )
                half_to_pair[barcode] = pair
            pair_to_cell[pair] = cell_id
            rows.append((cell_id, polydt, ranhex))
    if not rows:
        raise ValueError("primer pair table is empty")
    return rows


def prepare_paired_primer_glm_data(
    alevin_dir,
    salmon_ref,
    pair_file,
    *,
    ec_design="weighted",
    regularization_target="theta",
    min_eq=5,
    min_half_umis=500,
    probability_file=None,
    weight_caches=None,
):
    """Prepare paired primer observations with one latent row per cell.

    Each retained biological cell contributes two equally weighted responses,
    ``0.5 C_p / sum(C_p)``. The matching design is the vertical stack of
    ``0.5 A_polydT`` and ``0.5 A_ranhex``. Both halves therefore estimate one
    shared abundance vector while retaining primer-specific observation bias.
    Only complete pairs meeting ``min_half_umis`` in both halves are retained.
    """
    if ec_design not in {"binary", "weighted"}:
        raise ValueError("paired primer fitting requires binary or weighted design")
    alevin_dir = Path(alevin_dir)
    features, membership = load_alevin_structure(alevin_dir)
    barcodes, counts = load_alevin_counts(alevin_dir)
    features = np.asarray(features)
    barcodes = np.asarray(barcodes)
    if counts.shape[0] != len(barcodes):
        raise ValueError("count rows do not match alevin barcode rows")

    barcode_to_row = {}
    for index, barcode in enumerate(barcodes):
        if barcode in barcode_to_row:
            raise ValueError(f"duplicate alevin barcode: {barcode}")
        barcode_to_row[barcode] = index
    pairs = _read_primer_pairs(pair_file)
    group_by_row = np.full(len(barcodes), -1, dtype=np.int8)
    complete = []
    seen_half_barcodes = set()
    raw_totals = np.asarray(counts.sum(axis=1)).ravel()
    for cell_id, polydt, ranhex in pairs:
        if polydt in seen_half_barcodes or ranhex in seen_half_barcodes:
            raise ValueError("a half-cell barcode occurs in more than one primer pair")
        seen_half_barcodes.update((polydt, ranhex))
        poly_row = barcode_to_row.get(polydt)
        hex_row = barcode_to_row.get(ranhex)
        if poly_row is not None:
            group_by_row[poly_row] = 0
        if hex_row is not None:
            group_by_row[hex_row] = 1
        if (
            poly_row is not None and hex_row is not None
            and raw_totals[poly_row] >= float(min_half_umis)
            and raw_totals[hex_row] >= float(min_half_umis)
        ):
            complete.append((cell_id, poly_row, hex_row))
    if not complete:
        raise ValueError("no complete primer pairs meet min_half_umis")

    poly_rows = np.asarray([row[1] for row in complete], dtype=np.int64)
    hex_rows = np.asarray([row[2] for row in complete], dtype=np.int64)
    aggregate = np.asarray(
        counts[poly_rows].sum(axis=0) + counts[hex_rows].sum(axis=0)
    ).ravel()
    transcript_count = aggregate @ membership
    ec_keep = aggregate >= float(min_eq)
    feature_keep = np.asarray(transcript_count).ravel() > 0

    transcript_lengths = sc_utils.get_transcript_lengths(Path(salmon_ref))
    if ec_design == "binary":
        base = sc_utils._column_normalize(membership.astype(float))
        phi_designs = [base, base]
    else:
        probability_file = Path(probability_file or (
            alevin_dir / "gene_eqclass_probs.tsv.gz"
        ))
        if weight_caches is None:
            weight_caches = [
                alevin_dir / "gene_eqclass_fixed_weights_polydt.npz",
                alevin_dir / "gene_eqclass_fixed_weights_ranhex.npz",
            ]
        phi_designs = sc_utils.grouped_ec_probability_matrices(
            probability_file,
            membership,
            group_by_row,
            2,
            cache_files=weight_caches,
        )

    filtered_features = features[feature_keep]
    _, weights = sc_utils.get_feature_weights(filtered_features, transcript_lengths)
    designs = []
    for phi_design in phi_designs:
        filtered = phi_design[ec_keep, :][:, feature_keep]
        designs.append(sc_utils.parameterize_glm_design(
            filtered,
            weights,
            regularization_target,
            normalize_columns=True,
        ))
    compatibility = sp.vstack(
        (0.5 * designs[0], 0.5 * designs[1]), format="csr"
    )
    raw_poly_counts = counts[poly_rows][:, ec_keep].tocsr()
    raw_hex_counts = counts[hex_rows][:, ec_keep].tocsr()
    poly_counts, poly_totals = _row_normalize(raw_poly_counts)
    hex_counts, hex_totals = _row_normalize(raw_hex_counts)
    paired_counts = sp.hstack(
        (0.5 * poly_counts, 0.5 * hex_counts), format="csr"
    )
    pair_ids = np.asarray([row[0] for row in complete])
    return PreparedGLMData(
        paired_counts,
        compatibility,
        pair_ids,
        filtered_features,
        poly_totals + hex_totals,
        metadata={
            "paired_primers": True,
            "primer_names": ["polydT", "ranhex"],
            "min_half_umis": float(min_half_umis),
            "half_umi_totals": np.column_stack((poly_totals, hex_totals)),
            "source_rows": np.column_stack((poly_rows, hex_rows)),
            "annotated_pair_count": len(pairs),
            "retained_pair_count": len(complete),
        },
        cv_raw_counts=sp.hstack(
            (raw_poly_counts, raw_hex_counts), format="csr"
        ),
    )


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

    cell_umi_totals = np.asarray(counts.sum(axis=1)).ravel()
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
        filtered_counts,
        compatibility,
        barcodes,
        filtered_features,
        cell_umi_totals,
    )


def sample_cells_by_count(counts, n_cells, seed=0, *, min_count=1, totals=None):
    """Select cells meeting a count threshold, optionally subsampling them."""
    if float(min_count) < 0:
        raise ValueError("min_count must be nonnegative")
    if totals is None:
        totals = np.asarray(counts.sum(axis=1)).ravel()
    else:
        totals = np.asarray(totals).ravel()
        if len(totals) != counts.shape[0]:
            raise ValueError("totals must have one value per count-matrix row")
    eligible = np.flatnonzero(totals >= float(min_count))
    if n_cells is None or int(n_cells) <= 0 or int(n_cells) >= len(eligible):
        return eligible
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(eligible, size=int(n_cells), replace=False))


def sample_nonempty_cells(counts, n_cells, seed=0):
    """Select a reproducible random subset of nonempty cell rows."""
    return sample_cells_by_count(counts, n_cells, seed=seed, min_count=1)


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


def _standard_count_fold_pairs(counts, n_folds=3, seed=0):
    validation_folds = split_count_folds(counts, n_folds=n_folds, seed=seed)
    return [
        (
            sum(
                (
                    fold
                    for index, fold in enumerate(validation_folds)
                    if index != fold_index
                ),
                start=sp.csr_matrix(counts.shape),
            ).tocsr(),
            validation_counts,
        )
        for fold_index, validation_counts in enumerate(validation_folds)
    ]


def _normalize_paired_primer_counts(counts):
    """Give each primer half equal mass within every nonempty cell."""
    counts = counts.tocsr()
    if counts.shape[1] % 2:
        raise ValueError("paired-primer count matrices need two equal-width halves")
    half = counts.shape[1] // 2
    polydt, _ = _row_normalize(counts[:, :half])
    ranhex, _ = _row_normalize(counts[:, half:])
    return sp.hstack((0.5 * polydt, 0.5 * ranhex), format="csr")


def paired_primer_count_fold_pairs(raw_counts, n_folds=3, seed=0):
    """Split raw molecules, then equalize primer mass within each CV partition."""
    raw_counts = raw_counts.tocsr()
    raw_pairs = _standard_count_fold_pairs(
        raw_counts, n_folds=n_folds, seed=seed
    )
    return [
        (
            _normalize_paired_primer_counts(training),
            _normalize_paired_primer_counts(validation),
        )
        for training, validation in raw_pairs
    ]


def _resolve_count_fold_pairs(counts, n_folds, seed, fold_pairs):
    pairs = (
        _standard_count_fold_pairs(counts, n_folds=n_folds, seed=seed)
        if fold_pairs is None
        else list(fold_pairs)
    )
    if len(pairs) != int(n_folds):
        raise ValueError("fold_pairs must contain n_folds train/validation pairs")
    for training, validation in pairs:
        if training.shape != counts.shape or validation.shape != counts.shape:
            raise ValueError("fold-pair matrices must match the response shape")
    return pairs


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
                         batch_cells=4096, power_iter=10, seed=0,
                         data_backend="auto"):
    """Return lambda_max or a rank-one line-search tau reference."""
    data = (
        counts if isinstance(counts, glm_solvers.SparseGLM)
        else glm_solvers.SparseGLM(
            counts, compatibility, device=device, batch_cells=batch_cells,
            data_backend=data_backend,
        )
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
    data_backend="auto",
    power_iter=10,
    fit_kwargs=None,
    warm_start=True,
    min_profile_active_fraction=0.9,
    min_profile_relative_variance=1e-6,
    fold_pairs=None,
):
    """Tune a scale-free lambda fraction or tau multiplier by count folds."""
    fit_kwargs = dict(fit_kwargs or {})
    fold_pairs = _resolve_count_fold_pairs(
        counts, n_folds, seed, fold_pairs
    )
    rows = []
    scales = []
    for fold_index, (training_counts, validation_counts) in enumerate(fold_pairs):
        training = (
            glm_solvers.SparseGLM(
                training_counts,
                compatibility,
                device=device,
                batch_cells=batch_cells,
                data_backend=data_backend,
            )
            if method == "admm_factorized" else training_counts
        )
        scale = hyperparameter_scale(
            training,
            None if isinstance(training, glm_solvers.SparseGLM)
            else compatibility,
            method,
            device=device,
            batch_cells=batch_cells,
            power_iter=power_iter,
            seed=seed + fold_index,
            data_backend=data_backend,
        )
        scales.append(scale)
        validation = glm_solvers.SparseGLM(
            validation_counts,
            compatibility,
            device=device,
            batch_cells=batch_cells,
            data_backend=data_backend,
        )
        path_multipliers = sorted(
            (float(value) for value in multipliers),
            reverse=method == "admm_factorized",
        )
        warm_factors = None
        for multiplier in path_multipliers:
            kwargs = dict(fit_kwargs)
            if method == "admm_factorized":
                kwargs["regularization"] = float(multiplier) * scale
            elif method == "frank_wolfe_penalized":
                kwargs["tau"] = float(multiplier) * scale
            if warm_start and warm_factors is not None:
                kwargs["initial_factors"] = warm_factors
            result = glm_solvers.fit_glm(
                training,
                None if isinstance(training, glm_solvers.SparseGLM)
                else compatibility,
                method,
                device=device,
                batch_cells=batch_cells,
                **kwargs,
            )
            validation_loss = validation.loss_for_factors(
                result.left, result.right
            ) / validation.n_nonempty
            profile = glm_solvers.factor_profile_diagnostics(
                result, batch_cells=batch_cells
            )
            rows.append(
                {
                    "fold": fold_index,
                    "multiplier": float(multiplier),
                    "scale": scale,
                    "value": float(multiplier) * scale,
                    "validation_loss_per_cell": validation_loss,
                    "iterations": result.diagnostics.get("iterations"),
                    "converged": result.diagnostics.get("converged"),
                    "warm_started": result.diagnostics.get("warm_started", False),
                    "warm_start_rank": result.diagnostics.get("warm_start_rank", 0),
                    "data_backend": result.diagnostics.get("data_backend"),
                    "cells_per_second": result.diagnostics.get("cells_per_second"),
                    "mean_epoch_seconds": float(np.mean(
                        result.diagnostics.get("epoch_seconds", [np.nan])
                    )),
                    "peak_cuda_memory_bytes": result.diagnostics.get(
                        "peak_cuda_memory_bytes"
                    ),
                    **profile,
                }
            )
            if warm_start:
                warm_factors = (
                    result.left.detach(), result.right.detach()
                )
            del result
            gc.collect()
            if device == "cuda" or (
                device == "auto" and glm_solvers._torch().cuda.is_available()
            ):
                glm_solvers._torch().cuda.empty_cache()
    means = {}
    standard_errors = {}
    candidate_converged = {}
    candidate_nondegenerate = {}
    mean_profile_variance = {}
    minimum_profile_active_fraction = {}
    for multiplier in multipliers:
        candidate_rows = [
            row for row in rows if row["multiplier"] == float(multiplier)
        ]
        values = [row["validation_loss_per_cell"] for row in candidate_rows]
        means[float(multiplier)] = float(np.mean(values))
        standard_errors[float(multiplier)] = float(
            np.std(values, ddof=1) / np.sqrt(len(values))
            if len(values) > 1 else 0.0
        )
        candidate_converged[float(multiplier)] = bool(
            candidate_rows and all(row.get("converged") for row in candidate_rows)
        )
        candidate_nondegenerate[float(multiplier)] = bool(
            candidate_rows and all(
                row.get("normalized_profile_active_fraction", 0.0)
                >= float(min_profile_active_fraction)
                and row.get("normalized_profile_relative_variance", 0.0)
                > float(min_profile_relative_variance)
                for row in candidate_rows
            )
        )
        mean_profile_variance[float(multiplier)] = float(np.mean([
            row.get("normalized_profile_relative_variance", 0.0)
            for row in candidate_rows
        ]))
        minimum_profile_active_fraction[float(multiplier)] = float(min(
            row.get("normalized_profile_active_fraction", 0.0)
            for row in candidate_rows
        ))
    best_multiplier = min(means, key=means.get)
    best_on_boundary = _best_on_open_boundary(
        method, multipliers, best_multiplier
    )
    return {
        "method": method,
        "warm_start": bool(warm_start),
        "regularization_path": sorted(
            (float(value) for value in multipliers),
            reverse=method == "admm_factorized",
        ),
        "n_folds": int(n_folds),
        "multipliers": [float(value) for value in multipliers],
        "fold_scales": scales,
        "mean_validation_loss": means,
        "validation_standard_error": standard_errors,
        "candidate_converged": candidate_converged,
        "candidate_nondegenerate": candidate_nondegenerate,
        "mean_profile_relative_variance": mean_profile_variance,
        "minimum_profile_active_fraction": minimum_profile_active_fraction,
        "best_multiplier": best_multiplier,
        "best_on_boundary": best_on_boundary,
        "fold_results": rows,
    }


def cross_validate_factorized_rank(
    counts,
    compatibility,
    ranks,
    *,
    n_folds=3,
    seed=0,
    device="auto",
    batch_cells=4096,
    data_backend="auto",
    fit_kwargs=None,
    selection_rule="one_standard_error",
    require_converged=False,
    require_nondegenerate=False,
    min_profile_active_fraction=0.9,
    min_profile_relative_variance=1e-6,
    _state=None,
    progress_callback=None,
    _return_state=False,
    fold_pairs=None,
):
    """Select unregularized nonnegative factor rank by molecule-count folds."""
    if selection_rule not in {"minimum", "one_standard_error"}:
        raise ValueError(f"invalid rank selection rule: {selection_rule}")
    ranks = sorted(set(int(value) for value in ranks))
    if not ranks or ranks[0] < 1:
        raise ValueError("rank candidates must be positive")
    fit_kwargs = dict(fit_kwargs or {})
    fit_kwargs.pop("rank", None)
    fold_pairs = (
        _resolve_count_fold_pairs(counts, n_folds, seed, fold_pairs)
        if _state is None else _state["fold_pairs"]
    )
    rows = [] if _state is None else list(_state["rows"])
    continued_factors = [None] * len(fold_pairs)
    for fold_index, (training_counts, validation_counts) in enumerate(fold_pairs):
        completed = {
            row["rank"] for row in rows if row["fold"] == fold_index
        }
        pending_ranks = [rank for rank in ranks if rank not in completed]
        if not pending_ranks:
            continued_factors[fold_index] = (
                None if _state is None else _state["warm_factors"][fold_index]
            )
            continue
        training = glm_solvers.SparseGLM(
            training_counts,
            compatibility,
            device=device,
            batch_cells=batch_cells,
            data_backend=data_backend,
        )
        validation = glm_solvers.SparseGLM(
            validation_counts,
            compatibility,
            device=device,
            batch_cells=batch_cells,
            data_backend=data_backend,
        )
        warm_factors = (
            None if _state is None else _state["warm_factors"][fold_index]
        )
        if warm_factors is not None and pending_ranks[0] <= warm_factors[0].shape[1]:
            raise ValueError("incremental rank candidates must increase")
        for rank in pending_ranks:
            kwargs = dict(fit_kwargs, rank=rank)
            if warm_factors is not None:
                kwargs["initial_factors"] = warm_factors
            result = glm_solvers.fit_glm(
                training,
                None,
                "factorized",
                device=device,
                batch_cells=batch_cells,
                **kwargs,
            )
            validation_loss = validation.loss_for_factors(
                result.left, result.right
            ) / validation.n_nonempty
            profile = glm_solvers.factor_profile_diagnostics(
                result, batch_cells=batch_cells
            )
            row = {
                "fold": fold_index,
                "rank": rank,
                "validation_loss_per_cell": validation_loss,
                "iterations": result.diagnostics.get("iterations"),
                "converged": result.diagnostics.get("converged"),
                "warm_started": result.diagnostics.get(
                    "warm_started", False
                ),
                "warm_start_rank": result.diagnostics.get(
                    "warm_start_rank", 0
                ),
                "data_backend": result.diagnostics.get("data_backend"),
                "cells_per_second": result.diagnostics.get("cells_per_second"),
                "mean_epoch_seconds": float(np.mean(
                    result.diagnostics.get("epoch_seconds", [np.nan])
                )),
                "minibatch_iterations": result.diagnostics.get(
                    "minibatch_iterations"
                ),
                "polish_iterations": result.diagnostics.get(
                    "polish_iterations"
                ),
                "peak_cuda_memory_bytes": result.diagnostics.get(
                    "peak_cuda_memory_bytes"
                ),
                **profile,
            }
            rows.append(row)
            if progress_callback is not None:
                progress_callback(dict(row))
            warm_factors = (result.left.detach(), result.right.detach())
            del result
            gc.collect()
            if device == "cuda" or (
                device == "auto" and glm_solvers._torch().cuda.is_available()
            ):
                glm_solvers._torch().cuda.empty_cache()
        continued_factors[fold_index] = (
            warm_factors[0].detach().cpu(), warm_factors[1].detach().cpu()
        )
        del training, validation
        gc.collect()
        if device == "cuda" or (
            device == "auto" and glm_solvers._torch().cuda.is_available()
        ):
            glm_solvers._torch().cuda.empty_cache()

    means = {}
    standard_errors = {}
    candidate_converged = {}
    candidate_nondegenerate = {}
    mean_profile_variance = {}
    minimum_profile_active_fraction = {}
    for rank in ranks:
        candidate_rows = [row for row in rows if row["rank"] == rank]
        values = [row["validation_loss_per_cell"] for row in candidate_rows]
        means[rank] = float(np.mean(values))
        standard_errors[rank] = float(
            np.std(values, ddof=1) / np.sqrt(len(values))
            if len(values) > 1 else 0.0
        )
        candidate_converged[rank] = bool(
            candidate_rows and all(row.get("converged") for row in candidate_rows)
        )
        candidate_nondegenerate[rank] = bool(
            candidate_rows and all(
                row.get("normalized_profile_active_fraction", 0.0)
                >= float(min_profile_active_fraction)
                and row.get("normalized_profile_relative_variance", 0.0)
                > float(min_profile_relative_variance)
                for row in candidate_rows
            )
        )
        mean_profile_variance[rank] = float(np.mean([
            row.get("normalized_profile_relative_variance", 0.0)
            for row in candidate_rows
        ]))
        minimum_profile_active_fraction[rank] = float(min(
            row.get("normalized_profile_active_fraction", 0.0)
            for row in candidate_rows
        ))

    eligible = [
        rank for rank in ranks
        if not require_converged or candidate_converged[rank]
    ]
    convergence_met = bool(eligible)
    if require_nondegenerate:
        eligible = [rank for rank in eligible if candidate_nondegenerate[rank]]
    minimum = min(eligible, key=means.get) if eligible else None
    threshold = None
    best_rank = None
    if minimum is not None:
        threshold = means[minimum]
        if selection_rule == "one_standard_error":
            threshold += standard_errors[minimum]
            best_rank = min(rank for rank in eligible if means[rank] <= threshold)
        else:
            best_rank = minimum
    report = {
        "method": "factorized",
        "tuning_parameter": "rank",
        "regularization": None,
        "warm_start": True,
        "rank_path": ranks,
        "ranks": ranks,
        "n_folds": int(n_folds),
        "mean_validation_loss": means,
        "validation_standard_error": standard_errors,
        "candidate_converged": candidate_converged,
        "candidate_nondegenerate": candidate_nondegenerate,
        "mean_profile_relative_variance": mean_profile_variance,
        "minimum_profile_active_fraction": minimum_profile_active_fraction,
        "selection_rule": selection_rule,
        "require_converged": bool(require_converged),
        "convergence_requirement_met": convergence_met,
        "require_nondegenerate": bool(require_nondegenerate),
        "nondegeneracy_requirement_met": bool(eligible),
        "minimum_loss_rank": minimum,
        "selection_threshold": threshold,
        "best_rank": best_rank,
        "best_on_boundary": best_rank == ranks[-1],
        "fold_results": rows,
    }
    if not _return_state:
        return report
    return report, {
        "fold_pairs": fold_pairs,
        "rows": rows,
        "warm_factors": continued_factors,
    }


def cross_validate_factorized_rank_adaptive(
    counts,
    compatibility,
    ranks,
    *,
    max_rank=256,
    max_grid_expansions=2,
    **kwargs,
):
    """Expand a small-to-large rank path without replaying completed ranks."""
    candidates = sorted(set(int(value) for value in ranks))
    if not candidates or candidates[0] < 1:
        raise ValueError("rank candidates must be positive")
    if int(max_grid_expansions) < 0:
        raise ValueError("max_grid_expansions must be nonnegative")
    if int(max_rank) < candidates[-1]:
        raise ValueError("max_rank must not be below the initial rank grid")
    report, state = cross_validate_factorized_rank(
        counts, compatibility, candidates, _return_state=True, **kwargs
    )
    history = [{"round": 0, "added_rank": None, "best_rank": report["best_rank"]}]
    expansions = 0
    while (
        report["best_rank"] == candidates[-1]
        and candidates[-1] < int(max_rank)
        and expansions < int(max_grid_expansions)
    ):
        candidate = min(2 * candidates[-1], int(max_rank))
        if candidate in candidates:
            break
        candidates.append(candidate)
        report, state = cross_validate_factorized_rank(
            counts,
            compatibility,
            candidates,
            _state=state,
            _return_state=True,
            **kwargs,
        )
        expansions += 1
        history.append(
            {
                "round": expansions,
                "added_rank": candidate,
                "best_rank": report["best_rank"],
            }
        )
    report.update(
        grid_expansions=expansions,
        max_grid_expansions=int(max_grid_expansions),
        max_rank=int(max_rank),
        grid_exhausted=bool(
            report["best_rank"] == candidates[-1]
            and (
                candidates[-1] == int(max_rank)
                or expansions == int(max_grid_expansions)
            )
        ),
        grid_history=history,
    )
    return report


def _apply_selection_rule(report, method, selection_rule, require_converged,
                          require_nondegenerate=False,
                          profile_variance_retention=0.9):
    if selection_rule not in {
        "minimum", "one_standard_error", "one_se_variance_retention",
    }:
        raise ValueError(f"invalid CV selection rule: {selection_rule}")
    if not 0 < float(profile_variance_retention) <= 1:
        raise ValueError("profile_variance_retention must be in (0, 1]")
    candidates = [float(value) for value in report["multipliers"]]
    eligible = [
        value for value in candidates
        if not require_converged or report["candidate_converged"].get(value, False)
    ]
    convergence_met = bool(eligible)
    if require_nondegenerate:
        eligible = [
            value for value in eligible
            if report["candidate_nondegenerate"].get(value, False)
        ]
    report["selection_rule"] = selection_rule
    report["require_converged"] = bool(require_converged)
    report["convergence_requirement_met"] = convergence_met
    report["require_nondegenerate"] = bool(require_nondegenerate)
    report["nondegeneracy_requirement_met"] = bool(eligible)
    if not eligible:
        report.update(
            minimum_loss_multiplier=None,
            selection_threshold=None,
            best_multiplier=None,
            best_on_boundary=False,
        )
        return report

    minimum = min(eligible, key=report["mean_validation_loss"].get)
    threshold = report["mean_validation_loss"][minimum]
    profile_threshold = None
    if selection_rule in {"one_standard_error", "one_se_variance_retention"}:
        threshold += report["validation_standard_error"][minimum]
        within_one_se = [
            value for value in eligible
            if report["mean_validation_loss"][value] <= threshold
        ]
    if selection_rule == "one_standard_error":
        regularization_order = sorted(
            within_one_se, reverse=method == "admm_factorized"
        )
        selected = regularization_order[0]
    elif selection_rule == "one_se_variance_retention":
        maximum_profile_variance = max(
            report["mean_profile_relative_variance"][value]
            for value in within_one_se
        )
        profile_threshold = (
            float(profile_variance_retention) * maximum_profile_variance
        )
        structure_eligible = [
            value for value in within_one_se
            if report["mean_profile_relative_variance"][value]
            >= profile_threshold
        ]
        selected = sorted(
            structure_eligible, reverse=method == "admm_factorized"
        )[0]
    else:
        selected = minimum
    report.update(
        minimum_loss_multiplier=minimum,
        selection_threshold=threshold,
        profile_variance_retention=(
            float(profile_variance_retention)
            if selection_rule == "one_se_variance_retention" else None
        ),
        profile_variance_threshold=profile_threshold,
        best_multiplier=selected,
        best_on_boundary=_best_on_open_boundary(
            method, candidates, selected
        ),
    )
    return report


def cross_validate_glm_adaptive_grid(
    counts,
    compatibility,
    method,
    multipliers,
    *,
    max_grid_expansions=4,
    grid_expansion_factor=4.0,
    selection_rule="minimum",
    require_converged=False,
    require_nondegenerate=False,
    profile_variance_retention=0.9,
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
    _apply_selection_rule(
        report, method, selection_rule, require_converged,
        require_nondegenerate, profile_variance_retention,
    )
    history = [
        {
            "round": 0,
            "added_multiplier": None,
            "best_multiplier": report["best_multiplier"],
            "boundary_direction": (
                _open_boundary_direction(method, candidates, report["best_multiplier"])
                if report["best_multiplier"] is not None else None
            ),
        }
    ]
    expansions = 0
    while expansions < int(max_grid_expansions):
        if report["best_multiplier"] is None:
            break
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
        candidates.append(candidate)
        candidates.sort()
        # Replay the path so a new strongest endpoint is always followed by
        # strong-to-weak continuation, and every candidate has the same path
        # dependence after grid expansion.
        report = cross_validate_glm(
            counts, compatibility, method, candidates, **kwargs
        )
        _apply_selection_rule(
            report, method, selection_rule, require_converged,
            require_nondegenerate, profile_variance_retention,
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
