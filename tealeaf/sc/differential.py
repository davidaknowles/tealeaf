"""Subisoform covariance and differential-splicing utilities."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import re

import numpy as np
import scipy.linalg
import scipy.optimize
import scipy.sparse as sp
import scipy.stats


_ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]+)"')


@dataclass(frozen=True)
class SpliceBlock:
    """Alternative transcript paths between two constitutive anchors."""

    gene_id: str
    block_id: str
    chromosome: str
    strand: str
    left_anchor: tuple[int, int] | None
    right_anchor: tuple[int, int] | None
    transcripts: tuple[str, ...]
    path_index: tuple[int, ...]
    path_signatures: tuple[tuple[tuple[int, int], ...], ...]

    @property
    def n_paths(self):
        return len(self.path_signatures)


@dataclass
class CovarianceResult:
    """Covariance together with identifiability diagnostics."""

    covariance: np.ndarray
    information_rank: int
    parameter_dimension: int
    identifiable: bool
    null_projection: float
    minimum_positive_eigenvalue: float


@dataclass
class PathFit:
    """Local splice-path fit around a supplied transcript baseline."""

    delta: np.ndarray
    path_logratios: np.ndarray
    path_proportions: np.ndarray
    theta: np.ndarray
    covariance: CovarianceResult
    converged: bool
    iterations: int
    objective: float


def decode_hierarchical_theta(
    state,
    gene_intercept,
    gene_loadings,
    isoform_intercept,
    isoform_loadings,
    transcript_gene,
):
    """Decode normalized transcript abundance for a block of cell states."""
    state = np.asarray(state, dtype=np.float32)
    gene_intercept = np.asarray(gene_intercept, dtype=np.float32)
    gene_loadings = np.asarray(gene_loadings, dtype=np.float32)
    isoform_intercept = np.asarray(isoform_intercept, dtype=np.float32)
    isoform_loadings = np.asarray(isoform_loadings, dtype=np.float32)
    transcript_gene = np.asarray(transcript_gene, dtype=np.int64)
    gene_logits = np.clip(
        gene_intercept[:, None] + gene_loadings @ state.T, -20, 20
    )
    gene_logits -= gene_logits.max(axis=0, keepdims=True)
    gene_abundance = np.exp(gene_logits)
    gene_abundance /= gene_abundance.sum(axis=0, keepdims=True)

    order = np.argsort(transcript_gene, kind="stable")
    sorted_gene = transcript_gene[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_gene)) + 1]
    logits = isoform_intercept[:, None] + isoform_loadings @ state.T
    sorted_logits = logits[order]
    maxima = np.maximum.reduceat(sorted_logits, starts, axis=0)
    shifted = sorted_logits - np.repeat(
        maxima, np.diff(np.r_[starts, len(order)]), axis=0
    )
    unnormalized = np.exp(shifted)
    sums = np.add.reduceat(unnormalized, starts, axis=0)
    proportions = unnormalized / np.repeat(
        sums, np.diff(np.r_[starts, len(order)]), axis=0
    )
    sorted_theta = proportions * gene_abundance[sorted_gene]
    theta = np.empty_like(sorted_theta)
    theta[order] = sorted_theta
    return theta


def aggregate_hierarchical_theta(
    state,
    parameters,
    group_index,
    weights,
    *,
    batch_cells=128,
):
    """Decode cells and form UMI-weighted group transcript proportions."""
    state = np.asarray(state, dtype=np.float32)
    group_index = np.asarray(group_index, dtype=np.int64)
    weights = np.asarray(weights, dtype=np.float64)
    if len(state) != len(group_index) or len(state) != len(weights):
        raise ValueError("state, group_index, and weights must have equal lengths")
    valid = (group_index >= 0) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        raise ValueError("no cells have a valid group and positive weight")
    n_groups = int(group_index[valid].max()) + 1
    n_transcripts = len(parameters["isoform_intercept"])
    aggregate = np.zeros((n_groups, n_transcripts), dtype=np.float64)
    totals = np.zeros(n_groups, dtype=np.float64)
    for start in range(0, len(state), int(batch_cells)):
        stop = min(start + int(batch_cells), len(state))
        keep = valid[start:stop]
        if not keep.any():
            continue
        theta = decode_hierarchical_theta(
            state[start:stop][keep],
            parameters["gene_intercept"],
            parameters["gene_loadings"],
            parameters["isoform_intercept"],
            parameters["isoform_loadings"],
            parameters["transcript_gene"],
        )
        batch_groups = group_index[start:stop][keep]
        batch_weights = weights[start:stop][keep]
        for group in np.unique(batch_groups):
            selected = batch_groups == group
            aggregate[group] += theta[:, selected] @ batch_weights[selected]
            totals[group] += batch_weights[selected].sum()
    aggregate /= np.maximum(totals[:, None], np.finfo(float).tiny)
    return aggregate, totals


def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def _merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return tuple((start, end) for start, end in merged)


def read_gtf_exons(path, allowed_transcripts=None):
    """Read exon chains for selected transcripts from a GTF.

    Coordinates are returned as zero-based, half-open intervals.
    """
    allowed = None if allowed_transcripts is None else set(allowed_transcripts)
    genes = {}
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = dict(_ATTRIBUTE.findall(fields[8]))
            transcript = attributes.get("transcript_id")
            gene = attributes.get("gene_id")
            if transcript is None or gene is None:
                continue
            if allowed is not None and transcript not in allowed:
                continue
            record = genes.setdefault(
                gene,
                {
                    "chromosome": fields[0],
                    "strand": fields[6],
                    "transcripts": {},
                },
            )
            if (
                record["chromosome"] != fields[0]
                or record["strand"] != fields[6]
            ):
                raise ValueError(f"gene {gene} spans chromosomes or strands")
            record["transcripts"].setdefault(transcript, []).append(
                (int(fields[3]) - 1, int(fields[4]))
            )
    for record in genes.values():
        record["transcripts"] = {
            transcript: _merge_intervals(exons)
            for transcript, exons in record["transcripts"].items()
        }
    return genes


def _constitutive_intervals(exon_chains):
    boundaries = sorted({
        coordinate
        for chain in exon_chains
        for interval in chain
        for coordinate in interval
    })
    constitutive = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        if start == end:
            continue
        if all(
            any(
                exon_start <= start and exon_end >= end
                for exon_start, exon_end in chain
            )
            for chain in exon_chains
        ):
            constitutive.append((start, end))
    return _merge_intervals(constitutive)


def _clipped_signature(chain, start, end, strand):
    intervals = []
    for exon_start, exon_end in chain:
        clipped_start = max(start, exon_start)
        clipped_end = min(end, exon_end)
        if clipped_start < clipped_end:
            intervals.append((clipped_start, clipped_end))
    intervals = tuple(intervals)
    return intervals if strand == "+" else tuple(reversed(intervals))


def build_splice_blocks(gtf_path, allowed_transcripts=None):
    """Construct alternative paths between consecutive constitutive regions.

    Terminal pseudo-anchors are included, so alternative first and last
    regions are represented. Genes without a shared constitutive region get
    one whole-gene block.
    """
    genes = read_gtf_exons(gtf_path, allowed_transcripts)
    blocks = []
    for gene, record in sorted(genes.items()):
        transcript_items = sorted(record["transcripts"].items())
        if len(transcript_items) < 2:
            continue
        transcripts = tuple(item[0] for item in transcript_items)
        chains = tuple(item[1] for item in transcript_items)
        anchors = _constitutive_intervals(chains)
        gene_start = min(start for chain in chains for start, _ in chain)
        gene_end = max(end for chain in chains for _, end in chain)
        candidate_blocks = []
        if anchors:
            candidate_blocks.append((None, anchors[0], gene_start, anchors[0][0]))
            for left, right in zip(anchors[:-1], anchors[1:]):
                candidate_blocks.append((left, right, left[1], right[0]))
            candidate_blocks.append((anchors[-1], None, anchors[-1][1], gene_end))
        else:
            candidate_blocks.append((None, None, gene_start, gene_end))
        if record["strand"] == "-":
            candidate_blocks.reverse()
        block_number = 0
        for left, right, start, end in candidate_blocks:
            if start >= end:
                continue
            signatures = tuple(
                _clipped_signature(chain, start, end, record["strand"])
                for chain in chains
            )
            unique = sorted(set(signatures))
            if len(unique) < 2:
                continue
            signature_to_path = {
                signature: index for index, signature in enumerate(unique)
            }
            path_index = tuple(signature_to_path[value] for value in signatures)
            block_number += 1
            blocks.append(
                SpliceBlock(
                    gene_id=gene,
                    block_id=f"{gene}:B{block_number}",
                    chromosome=record["chromosome"],
                    strand=record["strand"],
                    left_anchor=left,
                    right_anchor=right,
                    transcripts=transcripts,
                    path_index=path_index,
                    path_signatures=tuple(unique),
                )
            )
    return blocks


def helmert_basis(size):
    """Return an orthonormal basis for the sum-zero subspace."""
    size = int(size)
    if size < 2:
        raise ValueError("a log-ratio basis requires at least two categories")
    basis = np.zeros((size, size - 1), dtype=float)
    for column in range(size - 1):
        scale = np.sqrt((column + 1) * (column + 2))
        basis[: column + 1, column] = 1.0 / scale
        basis[column + 1, column] = -(column + 1) / scale
    return basis


def transcript_fisher_information(theta, designs, totals):
    """Expected information for transcript log abundance.

    Each design is an EC-by-transcript matrix for one primer and each total is
    the number of gene-assigned molecules from that primer.
    """
    theta = np.asarray(theta, dtype=float)
    information = np.zeros((len(theta), len(theta)), dtype=float)
    for design, total in zip(designs, totals):
        total = float(total)
        if total <= 0:
            continue
        weighted = sp.csr_matrix(design).multiply(theta).tocsr()
        q = np.asarray(weighted.sum(axis=1)).ravel()
        positive = q > 0
        if not positive.any():
            continue
        weighted = weighted[positive]
        q = q[positive]
        normalizer = float(q.sum())
        if normalizer <= 0:
            continue
        cross = weighted.T @ sp.diags(1.0 / q) @ weighted
        cross = np.asarray(cross.toarray(), dtype=float) / normalizer
        exposure = np.asarray(weighted.sum(axis=0)).ravel() / normalizer
        information += total * (cross - np.outer(exposure, exposure))
    return 0.5 * (information + information.T)


def path_proportions(theta, path_index):
    """Collapse transcript abundance to normalized splice-path proportions."""
    theta = np.asarray(theta, dtype=float)
    path_index = np.asarray(path_index, dtype=np.int64)
    selected = path_index >= 0
    if not selected.any():
        raise ValueError("path_index does not select any transcripts")
    n_paths = int(path_index[selected].max()) + 1
    masses = np.bincount(
        path_index[selected], weights=theta[selected], minlength=n_paths
    ).astype(float)
    if np.any(masses <= 0):
        raise ValueError("every path needs positive baseline abundance")
    return masses / masses.sum()


def path_logratio_jacobian(theta, path_index, basis=None):
    """Jacobian of path ILR coordinates with respect to transcript log theta."""
    theta = np.asarray(theta, dtype=float)
    path_index = np.asarray(path_index, dtype=np.int64)
    selected = path_index >= 0
    n_paths = int(path_index[selected].max()) + 1
    basis = helmert_basis(n_paths) if basis is None else np.asarray(basis, dtype=float)
    proportions = path_proportions(theta, path_index)
    transcript_proportions = np.zeros_like(theta)
    transcript_proportions[selected] = theta[selected] / theta[selected].sum()
    derivative = np.zeros((n_paths, len(theta)), dtype=float)
    for transcript in np.flatnonzero(selected):
        path = path_index[transcript]
        derivative[:, transcript] = -transcript_proportions[transcript]
        derivative[path, transcript] += (
            transcript_proportions[transcript] / proportions[path]
        )
    return basis.T @ derivative


def identifiable_covariance(information, jacobian, rtol=1e-8, atol=1e-9):
    """Propagate inverse information while detecting unidentifiable contrasts."""
    information = np.asarray(information, dtype=float)
    jacobian = np.atleast_2d(np.asarray(jacobian, dtype=float))
    eigenvalues, eigenvectors = scipy.linalg.eigh(
        0.5 * (information + information.T), check_finite=False
    )
    scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
    positive = eigenvalues > scale * float(rtol)
    null = ~positive
    null_projection = (
        float(np.linalg.norm(jacobian @ eigenvectors[:, null]))
        if null.any()
        else 0.0
    )
    reference = max(float(np.linalg.norm(jacobian)), 1.0)
    identifiable = null_projection <= float(atol) * reference
    inverse = np.zeros_like(information)
    if positive.any():
        vectors = eigenvectors[:, positive]
        inverse = (vectors / eigenvalues[positive]) @ vectors.T
    covariance = jacobian @ inverse @ jacobian.T
    if not identifiable:
        covariance.fill(np.inf)
    return CovarianceResult(
        covariance=0.5 * (covariance + covariance.T),
        information_rank=int(positive.sum()),
        parameter_dimension=int(information.shape[0]),
        identifiable=bool(identifiable),
        null_projection=null_projection,
        minimum_positive_eigenvalue=(
            float(eigenvalues[positive].min()) if positive.any() else 0.0
        ),
    )


def profiled_path_covariance(theta, path_index, designs, totals, rtol=1e-8):
    """Path-logratio covariance after profiling all transcript directions."""
    information = transcript_fisher_information(theta, designs, totals)
    jacobian = path_logratio_jacobian(theta, path_index)
    return identifiable_covariance(information, jacobian, rtol=rtol)


def _perturbed_theta(baseline, path_index, basis, delta):
    baseline = np.asarray(baseline, dtype=float)
    path_index = np.asarray(path_index, dtype=np.int64)
    selected = path_index >= 0
    result = baseline.copy()
    mass = baseline[selected].sum()
    features = basis[path_index[selected]]
    logits = np.log(np.maximum(baseline[selected], np.finfo(float).tiny))
    logits += features @ delta
    logits -= logits.max()
    shares = np.exp(logits)
    shares /= shares.sum()
    result[selected] = mass * shares
    mean_feature = shares @ features
    log_jacobian = features - mean_feature
    return result, log_jacobian


def conditional_path_information(theta, path_index, basis, designs, totals):
    """Information for path perturbations with within-path mixtures fixed."""
    theta = np.asarray(theta, dtype=float)
    path_index = np.asarray(path_index, dtype=np.int64)
    selected = path_index >= 0
    features = basis[path_index[selected]]
    shares = theta[selected] / theta[selected].sum()
    log_jacobian = features - shares @ features
    dimension = basis.shape[1]
    information = np.zeros((dimension, dimension), dtype=float)
    for design, total in zip(designs, totals):
        total = float(total)
        if total <= 0:
            continue
        design = sp.csr_matrix(design)
        q = np.asarray(design @ theta).ravel()
        positive = q > 0
        if not positive.any():
            continue
        derivative = np.asarray(
            design[:, selected] @ (theta[selected, None] * log_jacobian)
        )
        q = q[positive]
        derivative = derivative[positive]
        normalizer = float(q.sum())
        normalizer_derivative = derivative.sum(axis=0) / normalizer
        information += total * (
            (derivative.T / q) @ derivative / normalizer
            - np.outer(normalizer_derivative, normalizer_derivative)
        )
    return 0.5 * (information + information.T)


def fit_path_perturbation(
    counts,
    designs,
    baseline,
    path_index,
    *,
    max_iter=100,
    tolerance=1e-7,
):
    """Fit local path logits and return conditional Fisher covariance."""
    baseline = np.maximum(np.asarray(baseline, dtype=float), 1e-12)
    path_index = np.asarray(path_index, dtype=np.int64)
    n_paths = int(path_index[path_index >= 0].max()) + 1
    basis = helmert_basis(n_paths)
    counts = tuple(np.asarray(value, dtype=float) for value in counts)
    designs = tuple(sp.csr_matrix(value) for value in designs)
    totals = tuple(float(value.sum()) for value in counts)
    baseline_paths = path_proportions(baseline, path_index)
    baseline_logratios = basis.T @ np.log(baseline_paths)

    def objective(delta):
        theta, log_jacobian = _perturbed_theta(
            baseline, path_index, basis, delta
        )
        selected = path_index >= 0
        loss = 0.0
        gradient = np.zeros_like(delta)
        for observed, design, total in zip(counts, designs, totals):
            if total <= 0:
                continue
            q = np.asarray(design @ theta).ravel()
            normalizer = float(q.sum())
            if normalizer <= 0 or np.any((observed > 0) & (q <= 0)):
                return np.inf, np.zeros_like(delta)
            loss -= float(observed @ np.log(np.maximum(q, 1e-300)))
            loss += total * np.log(normalizer)
            score = theta * np.asarray(
                design.T @ (observed / np.maximum(q, 1e-300))
            ).ravel()
            score -= total * theta * np.asarray(design.sum(axis=0)).ravel() / normalizer
            gradient -= log_jacobian.T @ score[selected]
        return loss, gradient

    initial_delta = np.zeros(n_paths - 1, dtype=float)
    if n_paths == 2:
        grid = np.linspace(-12.0, 12.0, 33)
        grid_losses = np.array([
            objective(np.array([value]))[0] for value in grid
        ])
        initial_delta[0] = grid[int(np.argmin(grid_losses))]
    result = scipy.optimize.minimize(
        objective,
        initial_delta,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": int(max_iter), "ftol": float(tolerance)},
    )
    theta, _ = _perturbed_theta(
        baseline, path_index, basis, np.asarray(result.x)
    )
    information = conditional_path_information(
        theta, path_index, basis, designs, totals
    )
    covariance = identifiable_covariance(
        information, np.eye(n_paths - 1), rtol=1e-8
    )
    proportions = path_proportions(theta, path_index)
    return PathFit(
        delta=np.asarray(result.x),
        path_logratios=basis.T @ np.log(proportions),
        path_proportions=proportions,
        theta=theta,
        covariance=covariance,
        converged=bool(result.success),
        iterations=int(result.nit),
        objective=float(result.fun),
    )


def multivariate_gls_test(
    values,
    covariances,
    design,
    tested_columns,
    *,
    biological_variance=None,
):
    """Precision-weighted Wald test with a REML scalar biological variance."""
    values = np.asarray(values, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    design = np.asarray(design, dtype=float)
    tested_columns = np.asarray(tested_columns, dtype=np.int64)
    if values.ndim != 2:
        raise ValueError("values must be samples by log-ratio dimensions")
    n_samples, dimension = values.shape
    if covariances.shape != (n_samples, dimension, dimension):
        raise ValueError("covariances have the wrong shape")
    if design.shape[0] != n_samples:
        raise ValueError("design and values have different sample counts")
    def fit_at_tau(tau2, model_design):
        parameter_count = model_design.shape[1] * dimension
        precision = np.zeros((parameter_count, parameter_count), dtype=float)
        score = np.zeros(parameter_count, dtype=float)
        log_determinant = 0.0
        quadratic_constant = 0.0
        for sample in range(n_samples):
            variance = covariances[sample] + tau2 * np.eye(dimension)
            sign, value = np.linalg.slogdet(variance)
            if sign <= 0:
                return None
            weight = scipy.linalg.inv(variance, check_finite=False)
            row_design = np.kron(
                model_design[sample : sample + 1], np.eye(dimension)
            )
            precision += row_design.T @ weight @ row_design
            score += row_design.T @ weight @ values[sample]
            quadratic_constant += values[sample] @ weight @ values[sample]
            log_determinant += value
        sign, precision_log_determinant = np.linalg.slogdet(precision)
        if sign <= 0:
            return None
        coefficient_covariance = scipy.linalg.inv(
            precision, check_finite=False
        )
        coefficients = coefficient_covariance @ score
        residual_quadratic = (
            quadratic_constant - score @ coefficient_covariance @ score
        )
        restricted_objective = (
            log_determinant
            + precision_log_determinant
            + max(float(residual_quadratic), 0.0)
        )
        return (
            restricted_objective,
            coefficients,
            coefficient_covariance,
        )

    null_design = np.delete(design, tested_columns, axis=1)
    if not null_design.shape[1]:
        raise ValueError("at least one untested design column is required")
    if biological_variance is None:
        initial = np.linalg.lstsq(null_design, values, rcond=None)[0]
        residual_scale = float(
            np.mean(np.square(values - null_design @ initial))
        )
        covariance_scale = float(
            np.mean(np.trace(covariances, axis1=1, axis2=2)) / dimension
        )
        scale = max(residual_scale, covariance_scale, 1e-8)

        def log_scale_objective(log_tau2):
            fitted = fit_at_tau(np.exp(log_tau2), null_design)
            return np.inf if fitted is None else fitted[0]

        optimized = scipy.optimize.minimize_scalar(
            log_scale_objective,
            bounds=(np.log(scale) - 18.0, np.log(scale) + 7.0),
            method="bounded",
            options={"xatol": 1e-5},
        )
        candidates = [(0.0, fit_at_tau(0.0, null_design))]
        if optimized.success:
            optimized_tau2 = float(np.exp(optimized.x))
            candidates.append(
                (
                    optimized_tau2,
                    fit_at_tau(optimized_tau2, null_design),
                )
            )
        candidates = [
            candidate for candidate in candidates if candidate[1] is not None
        ]
        if not candidates:
            raise np.linalg.LinAlgError("GLS null design has singular precision")
        tau2 = min(candidates, key=lambda candidate: candidate[1][0])[0]
    else:
        tau2 = max(float(biological_variance), 0.0)
    fitted = fit_at_tau(tau2, design)
    if fitted is None:
        raise np.linalg.LinAlgError("GLS design has singular precision")
    _, coefficients, coefficient_covariance = fitted
    tested = np.concatenate([
        np.arange(column * dimension, (column + 1) * dimension)
        for column in tested_columns
    ])
    tested_values = coefficients[tested]
    tested_covariance = coefficient_covariance[np.ix_(tested, tested)]
    statistic = float(
        tested_values
        @ scipy.linalg.pinvh(tested_covariance, rtol=1e-10)
        @ tested_values
    )
    degrees_of_freedom = int(np.linalg.matrix_rank(tested_covariance))
    return {
        "statistic": statistic,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": float(scipy.stats.chi2.sf(statistic, degrees_of_freedom)),
        "biological_variance": tau2,
        "coefficients": coefficients.reshape(design.shape[1], dimension),
        "coefficient_covariance": coefficient_covariance,
    }
