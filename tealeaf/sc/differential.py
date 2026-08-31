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
import scipy.special
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


def _block_value(block, name):
    if isinstance(block, dict):
        return block[name]
    return getattr(block, name)


def _path_pair_event(first, second, strand):
    """Classify the structural difference between two local exon chains."""
    first = tuple(tuple(interval) for interval in first)
    second = tuple(tuple(interval) for interval in second)
    first_unique = tuple(interval for interval in first if interval not in second)
    second_unique = tuple(interval for interval in second if interval not in first)
    if not first_unique or not second_unique:
        inserted = second_unique if not first_unique else first_unique
        if len(inserted) == 1:
            return "cassette exon"
        if len(inserted) > 1:
            return "multi-exon skipping"
    if len(first_unique) == len(second_unique) == 1:
        (first_start, first_end), (second_start, second_end) = first_unique[0], second_unique[0]
        if first_start == second_start and first_end != second_end:
            return "alternative donor" if strand == "+" else "alternative acceptor"
        if first_end == second_end and first_start != second_start:
            return "alternative acceptor" if strand == "+" else "alternative donor"
        if first_end <= second_start or second_end <= first_start:
            return "mutually exclusive exons"
    for one, many in ((first_unique, second_unique), (second_unique, first_unique)):
        many_start = min((interval[0] for interval in many), default=0)
        many_end = max((interval[1] for interval in many), default=0)
        if len(one) == 1 and len(many) > 1 and one[0][0] <= many_start and one[0][1] >= many_end:
            return "retained intron"
    return "complex internal"


def classify_splice_block(block, path_signatures=None):
    """Return a mutually exclusive structural event class for a splice block.

    ``path_signatures`` can restrict annotation paths to the paths represented
    by a fitted test. Terminal classes are called regions because one terminal
    block can contain several alternative exons or transcript ends.
    """
    signatures = _block_value(block, "path_signatures") if path_signatures is None else path_signatures
    signatures = tuple(tuple(tuple(interval) for interval in path) for path in signatures)
    left_anchor = _block_value(block, "left_anchor")
    right_anchor = _block_value(block, "right_anchor")
    strand = _block_value(block, "strand")
    if left_anchor is None and right_anchor is None:
        return "whole-gene complex"
    if left_anchor is None:
        return "alternative first region" if strand == "+" else "alternative last region"
    if right_anchor is None:
        return "alternative last region" if strand == "+" else "alternative first region"
    if len(signatures) == 2:
        return _path_pair_event(signatures[0], signatures[1], strand)
    return "compound internal"


def splice_block_event_components(block, path_signatures=None):
    """Return canonical event components found among a block's path pairs."""
    signatures = _block_value(block, "path_signatures") if path_signatures is None else path_signatures
    signatures = tuple(tuple(tuple(interval) for interval in path) for path in signatures)
    strand = _block_value(block, "strand")
    components = {
        _path_pair_event(signatures[first], signatures[second], strand)
        for first in range(len(signatures))
        for second in range(first + 1, len(signatures))
    }
    components.discard("complex internal")
    return tuple(sorted(components))


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


def path_laplace_log_evidence(
    fit, path_pseudocount, *, prior_center=None
):
    """Return the Laplace evidence for one ILR smoothing prior.

    The prior density with respect to orthonormal ILR coordinates is
    proportional to ``prod(psi**alpha)``. Constants independent of ``alpha``
    are omitted, which is sufficient for empirical-Bayes selection within a
    fixed path dimension.
    """
    alpha = float(path_pseudocount)
    if alpha <= 0:
        raise ValueError("path_pseudocount must be positive for evidence")
    n_paths = len(fit.path_proportions)
    dimension = n_paths - 1
    if not fit.covariance.identifiable:
        return -np.inf
    sign, covariance_log_determinant = np.linalg.slogdet(
        fit.covariance.covariance
    )
    if sign <= 0:
        return -np.inf
    if prior_center is None:
        prior_center = np.full(n_paths, 1.0 / n_paths)
    prior_center = np.asarray(prior_center, dtype=float)
    if prior_center.shape != (n_paths,) or np.any(prior_center <= 0):
        raise ValueError("prior_center must be a positive path composition")
    prior_center = prior_center / prior_center.sum()
    prior_counts = alpha * n_paths * prior_center
    log_prior_normalizer = (
        scipy.special.gammaln(prior_counts.sum())
        - scipy.special.gammaln(prior_counts).sum()
    )
    return float(
        -fit.objective
        + log_prior_normalizer
        + 0.5 * dimension * np.log(2 * np.pi)
        + 0.5 * covariance_log_determinant
    )


@dataclass
class SharedPathFit:
    """Shared path shift fit to cells with distinct transcript baselines."""

    delta: np.ndarray
    path_logratios: np.ndarray
    path_proportions: np.ndarray
    theta: np.ndarray
    covariance: CovarianceResult
    parameter_covariance: CovarianceResult
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


def _perturbed_theta_matrix(baseline, path_index, basis, delta):
    baseline = np.maximum(np.asarray(baseline, dtype=float), 1e-12)
    path_index = np.asarray(path_index, dtype=np.int64)
    selected = path_index >= 0
    result = baseline.copy()
    masses = baseline[:, selected].sum(axis=1)
    features = basis[path_index[selected]]
    logits = np.log(baseline[:, selected]) + features @ delta
    logits -= logits.max(axis=1, keepdims=True)
    shares = np.exp(logits)
    shares /= shares.sum(axis=1, keepdims=True)
    result[:, selected] = masses[:, None] * shares
    mean_features = shares @ features
    log_jacobian = features[None, :, :] - mean_features[:, None, :]
    return result, log_jacobian


def _shared_path_information(
    theta,
    log_jacobian,
    path_index,
    designs,
    counts,
):
    selected = np.asarray(path_index, dtype=np.int64) >= 0
    dimension = log_jacobian.shape[2]
    information = np.zeros((dimension, dimension), dtype=float)
    for design, observed in zip(designs, counts):
        design = sp.csr_matrix(design)
        observed = sp.csr_matrix(observed)
        totals = np.asarray(observed.sum(axis=1)).ravel()
        q = np.asarray(design @ theta.T).T
        normalizers = q.sum(axis=1)
        derivatives = np.empty(
            (len(theta), design.shape[0], dimension),
            dtype=float,
        )
        for coordinate in range(dimension):
            transcript_derivative = (
                theta[:, selected] * log_jacobian[:, :, coordinate]
            )
            derivatives[:, :, coordinate] = np.asarray(
                design[:, selected] @ transcript_derivative.T
            ).T
        positive = (totals > 0) & (normalizers > 0)
        if not positive.any():
            continue
        q = np.maximum(q[positive], 1e-300)
        derivatives = derivatives[positive]
        totals = totals[positive]
        normalizers = normalizers[positive]
        weights = totals[:, None] / (normalizers[:, None] * q)
        information += np.einsum(
            "ced,cef,ce->df",
            derivatives,
            derivatives,
            weights,
        )
        normalizer_gradient = (
            derivatives.sum(axis=1) / normalizers[:, None]
        )
        information -= np.einsum(
            "c,cd,ce->de",
            totals,
            normalizer_gradient,
            normalizer_gradient,
        )
    return 0.5 * (information + information.T)


def _weighted_path_response(
    theta,
    log_jacobian,
    path_index,
    basis,
    weights,
):
    path_index = np.asarray(path_index, dtype=np.int64)
    weights = np.asarray(weights, dtype=float)
    selected = path_index >= 0
    n_paths = basis.shape[0]
    masses = np.zeros(n_paths, dtype=float)
    derivatives = np.zeros((n_paths, basis.shape[1]), dtype=float)
    for local_transcript, transcript in enumerate(np.flatnonzero(selected)):
        path = path_index[transcript]
        weighted_theta = weights * theta[:, transcript]
        masses[path] += weighted_theta.sum()
        derivatives[path] += (
            weighted_theta[:, None]
            * log_jacobian[:, local_transcript, :]
        ).sum(axis=0)
    if np.any(masses <= 0):
        raise ValueError("every path needs positive weighted abundance")
    proportions = masses / masses.sum()
    response = basis.T @ np.log(proportions)
    jacobian = basis.T @ (derivatives / masses[:, None])
    return proportions, response, jacobian


def fit_shared_path_perturbation(
    counts,
    designs,
    baselines,
    path_index,
    weights,
    *,
    max_iter=100,
    tolerance=1e-7,
):
    """Fit one path shift across cells with cell-specific transcript baselines."""
    baselines = np.maximum(np.asarray(baselines, dtype=float), 1e-12)
    path_index = np.asarray(path_index, dtype=np.int64)
    weights = np.asarray(weights, dtype=float)
    if baselines.ndim != 2 or len(baselines) != len(weights):
        raise ValueError("baselines must be cells by transcripts")
    n_paths = int(path_index[path_index >= 0].max()) + 1
    basis = helmert_basis(n_paths)
    counts = tuple(sp.csr_matrix(value) for value in counts)
    designs = tuple(sp.csr_matrix(value) for value in designs)

    def objective(delta):
        theta, log_jacobian = _perturbed_theta_matrix(
            baselines, path_index, basis, delta
        )
        selected = path_index >= 0
        loss = 0.0
        gradient = np.zeros_like(delta)
        for observed, design in zip(counts, designs):
            totals = np.asarray(observed.sum(axis=1)).ravel()
            q = np.asarray(design @ theta.T).T
            normalizers = q.sum(axis=1)
            nonzero = observed.tocoo()
            positive_counts = nonzero.data > 0
            rows = nonzero.row[positive_counts]
            columns = nonzero.col[positive_counts]
            observed_values = nonzero.data[positive_counts]
            selected_q = q[rows, columns]
            if np.any(selected_q <= 0) or np.any(
                (totals > 0) & (normalizers <= 0)
            ):
                return np.inf, np.zeros_like(delta)
            loss -= float(observed_values @ np.log(selected_q))
            loss += float(
                totals @ np.log(np.maximum(normalizers, 1e-300))
            )
            ratio = observed.multiply(1.0 / np.maximum(q, 1e-300))
            score = np.asarray((ratio @ design).toarray()) * theta
            column_sums = np.asarray(design.sum(axis=0)).ravel()
            score -= (
                totals / np.maximum(normalizers, 1e-300)
            )[:, None] * theta * column_sums[None, :]
            gradient -= np.einsum(
                "cs,csd->d",
                score[:, selected],
                log_jacobian,
            )
        return loss, gradient

    initial_delta = np.zeros(n_paths - 1, dtype=float)
    if n_paths == 2:
        grid = np.linspace(-12.0, 12.0, 33)
        losses = np.array([
            objective(np.array([value]))[0] for value in grid
        ])
        initial_delta[0] = grid[int(np.argmin(losses))]
    result = scipy.optimize.minimize(
        objective,
        initial_delta,
        method="L-BFGS-B",
        jac=True,
        bounds=[(-20.0, 20.0)] * (n_paths - 1),
        options={"maxiter": int(max_iter), "ftol": float(tolerance)},
    )
    theta, log_jacobian = _perturbed_theta_matrix(
        baselines,
        path_index,
        basis,
        np.asarray(result.x),
    )
    information = _shared_path_information(
        theta,
        log_jacobian,
        path_index,
        designs,
        counts,
    )
    parameter_covariance = identifiable_covariance(
        information,
        np.eye(n_paths - 1),
    )
    proportions, response, response_jacobian = _weighted_path_response(
        theta,
        log_jacobian,
        path_index,
        basis,
        weights,
    )
    covariance = identifiable_covariance(
        information,
        response_jacobian,
    )
    return SharedPathFit(
        delta=np.asarray(result.x),
        path_logratios=response,
        path_proportions=proportions,
        theta=theta,
        covariance=covariance,
        parameter_covariance=parameter_covariance,
        converged=bool(result.success),
        iterations=int(result.nit),
        objective=float(result.fun),
    )


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
    path_pseudocount=0.0,
    path_prior_center="uniform",
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
    path_pseudocount = float(path_pseudocount)
    if path_pseudocount < 0:
        raise ValueError("path_pseudocount must be nonnegative")
    if path_prior_center == "uniform":
        prior_center = np.full(n_paths, 1.0 / n_paths)
    elif path_prior_center == "baseline":
        prior_center = baseline_paths
    else:
        raise ValueError("path_prior_center must be uniform or baseline")
    prior_counts = path_pseudocount * n_paths * prior_center

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
        if path_pseudocount:
            proportions = path_proportions(theta, path_index)
            loss -= float(
                prior_counts @ np.log(np.maximum(proportions, 1e-300))
            )
            gradient += (
                prior_counts.sum() * (proportions @ basis)
                - prior_counts @ basis
            )
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
        bounds=[(-20.0, 20.0)] * (n_paths - 1),
        options={"maxiter": int(max_iter), "ftol": float(tolerance)},
    )
    theta, _ = _perturbed_theta(
        baseline, path_index, basis, np.asarray(result.x)
    )
    information = conditional_path_information(
        theta, path_index, basis, designs, totals
    )
    if path_pseudocount:
        proportions = path_proportions(theta, path_index)
        simplex_covariance = (
            np.diag(proportions) - np.outer(proportions, proportions)
        )
        information += (
            prior_counts.sum()
            * basis.T
            @ simplex_covariance
            @ basis
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


def paired_measurement_error_test(
    differences,
    covariances,
    *,
    uncertainty_scale=1.0,
    biological_variance=None,
):
    """Test a paired mean while retaining conditional measurement covariance.

    ``differences`` has shape ``(M, D)`` and ``covariances`` has shape
    ``(M, D, D)``, where ``M`` is the number of paired subjects and ``D`` is
    the number of path log-ratio coordinates. A scalar isotropic biological
    variance is estimated by restricted maximum likelihood. The returned
    Wald reference is intended to be calibrated by design-matched sign flips.
    """
    differences = np.asarray(differences, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    if differences.ndim != 2 or differences.shape[1] < 1:
        raise ValueError("differences must be subjects by coordinates")
    n_subjects, dimension = differences.shape
    if covariances.shape != (n_subjects, dimension, dimension):
        raise ValueError("covariances have the wrong shape")
    uncertainty_scale = float(uncertainty_scale)
    if not np.isfinite(uncertainty_scale) or uncertainty_scale < 0:
        raise ValueError("uncertainty_scale must be finite and nonnegative")
    covariances = uncertainty_scale * covariances
    if n_subjects < 3 or n_subjects <= dimension:
        return {
            "statistic": 0.0,
            "degrees_of_freedom": dimension,
            "p_value": 1.0,
            "n_subjects": n_subjects,
            "biological_variance": np.nan,
            "mean": np.zeros(dimension),
            "mean_covariance": np.full((dimension, dimension), np.inf),
            "converged": False,
        }

    identity = np.eye(dimension)

    def fit_at_tau(tau2):
        precision = np.zeros((dimension, dimension), dtype=float)
        score = np.zeros(dimension, dtype=float)
        log_determinant = 0.0
        quadratic_constant = 0.0
        for value, covariance in zip(differences, covariances):
            variance = covariance + float(tau2) * identity
            sign, local_log_determinant = np.linalg.slogdet(variance)
            if sign <= 0:
                return None
            weight = scipy.linalg.inv(variance, check_finite=False)
            precision += weight
            score += weight @ value
            quadratic_constant += value @ weight @ value
            log_determinant += local_log_determinant
        sign, precision_log_determinant = np.linalg.slogdet(precision)
        if sign <= 0:
            return None
        mean_covariance = scipy.linalg.inv(precision, check_finite=False)
        mean = mean_covariance @ score
        residual_quadratic = quadratic_constant - score @ mean
        restricted_objective = (
            log_determinant
            + precision_log_determinant
            + max(float(residual_quadratic), 0.0)
        )
        return restricted_objective, mean, mean_covariance

    empirical_scale = max(
        float(np.mean(np.var(differences, axis=0, ddof=1))),
        float(np.mean(np.trace(covariances, axis1=1, axis2=2)) / dimension),
        1e-8,
    )
    if biological_variance is None:
        optimized = scipy.optimize.minimize_scalar(
            lambda log_tau2: (
                np.inf
                if (fit := fit_at_tau(np.exp(log_tau2))) is None
                else fit[0]
            ),
            bounds=(np.log(empirical_scale) - 18.0, np.log(empirical_scale) + 7.0),
            method="bounded",
            options={"xatol": 1e-5},
        )
        candidates = [(0.0, fit_at_tau(0.0))]
        if optimized.success:
            tau2 = float(np.exp(optimized.x))
            candidates.append((tau2, fit_at_tau(tau2)))
        candidates = [item for item in candidates if item[1] is not None]
    else:
        tau2 = max(float(biological_variance), 0.0)
        candidates = [(tau2, fit_at_tau(tau2))]
    if not candidates:
        raise np.linalg.LinAlgError("paired measurement covariance is singular")
    biological_variance, fitted = min(
        candidates, key=lambda item: item[1][0]
    )
    restricted_objective, mean, mean_covariance = fitted
    statistic = float(
        mean
        @ scipy.linalg.pinvh(mean_covariance, rtol=1e-10)
        @ mean
    )
    return {
        "statistic": statistic,
        "degrees_of_freedom": dimension,
        "p_value": float(scipy.stats.chi2.sf(statistic, dimension)),
        "n_subjects": n_subjects,
        "biological_variance": biological_variance,
        "restricted_objective": restricted_objective,
        "uncertainty_scale": uncertainty_scale,
        "mean": mean,
        "mean_covariance": mean_covariance,
        "converged": True,
    }


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
    if np.count_nonzero(covariances) == 0:
        null_design = np.delete(design, tested_columns, axis=1)
        if not null_design.shape[1]:
            raise ValueError("at least one untested design column is required")
        null_coefficients = np.linalg.lstsq(
            null_design, values, rcond=None
        )[0]
        null_residuals = values - null_design @ null_coefficients
        residual_degrees = (
            n_samples - np.linalg.matrix_rank(null_design)
        ) * dimension
        if residual_degrees <= 0:
            raise ValueError("GLS null model has no residual degrees of freedom")
        if biological_variance is None:
            residual_scale = float(np.mean(np.square(null_residuals)))
            lower_variance = max(residual_scale, 1e-8) * np.exp(-18.0)
            tau2 = max(
                float(np.sum(np.square(null_residuals)) / residual_degrees),
                lower_variance,
            )
        else:
            tau2 = max(float(biological_variance), np.finfo(float).tiny)
        null_crossproduct = null_design.T @ null_design
        sign, null_log_determinant = np.linalg.slogdet(null_crossproduct)
        if sign <= 0:
            raise np.linalg.LinAlgError("GLS null design has singular precision")
        null_restricted_objective = float(
            residual_degrees * np.log(tau2)
            + dimension * null_log_determinant
            + np.sum(np.square(null_residuals)) / tau2
        )
        crossproduct = design.T @ design
        sign, _ = np.linalg.slogdet(crossproduct)
        if sign <= 0:
            raise np.linalg.LinAlgError("GLS design has singular precision")
        crossproduct_inverse = scipy.linalg.inv(
            crossproduct, check_finite=False
        )
        coefficients = crossproduct_inverse @ design.T @ values
        coefficient_covariance = tau2 * np.kron(
            crossproduct_inverse, np.eye(dimension)
        )
        tested = np.concatenate([
            np.arange(column * dimension, (column + 1) * dimension)
            for column in tested_columns
        ])
        tested_values = coefficients[tested_columns].ravel()
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
            "p_value": float(
                scipy.stats.chi2.sf(statistic, degrees_of_freedom)
            ),
            "biological_variance": tau2,
            "restricted_objective": null_restricted_objective,
            "coefficients": coefficients,
            "coefficient_covariance": coefficient_covariance,
        }
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
        selected = min(candidates, key=lambda candidate: candidate[1][0])
        tau2 = selected[0]
        null_restricted_objective = float(selected[1][0])
    else:
        tau2 = max(float(biological_variance), 0.0)
        null_fit = fit_at_tau(tau2, null_design)
        if null_fit is None:
            raise np.linalg.LinAlgError("GLS null design has singular precision")
        null_restricted_objective = float(null_fit[0])
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
        "restricted_objective": null_restricted_objective,
        "coefficients": coefficients.reshape(design.shape[1], dimension),
        "coefficient_covariance": coefficient_covariance,
    }


def clustered_multivariate_gls_test(
    values,
    covariances,
    design,
    tested_columns,
    clusters,
    *,
    variance_components=None,
):
    """GLS test with shared-cluster and observation variance components."""
    values = np.asarray(values, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    design = np.asarray(design, dtype=float)
    tested_columns = np.asarray(tested_columns, dtype=np.int64)
    clusters = np.asarray(clusters)
    if values.ndim != 2:
        raise ValueError("values must be samples by log-ratio dimensions")
    n_samples, dimension = values.shape
    if covariances.shape != (n_samples, dimension, dimension):
        raise ValueError("covariances have the wrong shape")
    if design.shape[0] != n_samples or len(clusters) != n_samples:
        raise ValueError("design, clusters, and values must align")
    cluster_rows = [
        np.flatnonzero(clusters == cluster)
        for cluster in np.unique(clusters)
    ]

    def fit_at_variances(cluster_variance, residual_variance, model_design):
        parameter_count = model_design.shape[1] * dimension
        precision = np.zeros((parameter_count, parameter_count), dtype=float)
        score = np.zeros(parameter_count, dtype=float)
        log_determinant = 0.0
        quadratic_constant = 0.0
        for rows in cluster_rows:
            size = len(rows)
            variance = scipy.linalg.block_diag(*[
                covariances[row] + residual_variance * np.eye(dimension)
                for row in rows
            ])
            variance += cluster_variance * np.kron(
                np.ones((size, size)), np.eye(dimension)
            )
            sign, value = np.linalg.slogdet(variance)
            if sign <= 0:
                return None
            weight = scipy.linalg.inv(variance, check_finite=False)
            block_design = np.kron(model_design[rows], np.eye(dimension))
            block_values = values[rows].ravel()
            precision += block_design.T @ weight @ block_design
            score += block_design.T @ weight @ block_values
            quadratic_constant += block_values @ weight @ block_values
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
        objective = (
            log_determinant
            + precision_log_determinant
            + max(float(residual_quadratic), 0.0)
        )
        return objective, coefficients, coefficient_covariance

    null_design = np.delete(design, tested_columns, axis=1)
    if not null_design.shape[1]:
        raise ValueError("at least one untested design column is required")

    def finalize(cluster_variance, residual_variance):
        fitted = fit_at_variances(
            cluster_variance,
            residual_variance,
            design,
        )
        if fitted is None:
            raise np.linalg.LinAlgError(
                "clustered GLS alternative fit failed"
            )
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
            "p_value": float(
                scipy.stats.chi2.sf(statistic, degrees_of_freedom)
            ),
            "cluster_variance": cluster_variance,
            "residual_variance": residual_variance,
            "coefficients": coefficients.reshape(
                design.shape[1], dimension
            ),
            "coefficient_covariance": coefficient_covariance,
        }

    if variance_components is not None:
        variance_components = np.asarray(variance_components, dtype=float)
        if (
            variance_components.shape != (2,)
            or not np.isfinite(variance_components).all()
            or np.any(variance_components < 0)
        ):
            raise ValueError(
                "variance_components must contain two nonnegative values"
            )
        return finalize(*variance_components)

    initial = np.linalg.lstsq(null_design, values, rcond=None)[0]
    residual_scale = float(
        np.mean(np.square(values - null_design @ initial))
    )
    covariance_scale = float(
        np.mean(np.trace(covariances, axis1=1, axis2=2)) / dimension
    )
    scale = max(residual_scale, covariance_scale, 1e-8)
    lower = np.log(scale) - 18.0
    upper = np.log(scale) + 7.0
    candidates = [
        (
            0.0,
            0.0,
            fit_at_variances(0.0, 0.0, null_design),
        )
    ]

    def axis_fit(cluster_axis):
        def objective(log_variance):
            variance = float(np.exp(log_variance))
            fitted = fit_at_variances(
                variance if cluster_axis else 0.0,
                0.0 if cluster_axis else variance,
                null_design,
            )
            return np.inf if fitted is None else fitted[0]

        result = scipy.optimize.minimize_scalar(
            objective,
            bounds=(lower, upper),
            method="bounded",
            options={"xatol": 1e-5},
        )
        if result.success:
            variance = float(np.exp(result.x))
            candidates.append((
                variance if cluster_axis else 0.0,
                0.0 if cluster_axis else variance,
                fit_at_variances(
                    variance if cluster_axis else 0.0,
                    0.0 if cluster_axis else variance,
                    null_design,
                ),
            ))

    axis_fit(True)
    axis_fit(False)

    def joint_objective(log_variances):
        fitted = fit_at_variances(
            np.exp(log_variances[0]),
            np.exp(log_variances[1]),
            null_design,
        )
        return np.inf if fitted is None else fitted[0]

    joint = scipy.optimize.minimize(
        joint_objective,
        np.log([scale / 2.0, scale / 2.0]),
        method="L-BFGS-B",
        bounds=[(lower, upper), (lower, upper)],
        options={"maxiter": 100, "ftol": 1e-9},
    )
    if joint.success:
        cluster_variance, residual_variance = np.exp(joint.x)
        candidates.append((
            float(cluster_variance),
            float(residual_variance),
            fit_at_variances(
                cluster_variance,
                residual_variance,
                null_design,
            ),
        ))
    candidates = [
        candidate for candidate in candidates if candidate[2] is not None
    ]
    if not candidates:
        raise np.linalg.LinAlgError("clustered GLS null fit failed")
    cluster_variance, residual_variance, _ = min(
        candidates, key=lambda candidate: candidate[2][0]
    )
    return finalize(cluster_variance, residual_variance)


def _reference_multinomial_mean(design, coefficients):
    """Return reference-category multinomial means and derivatives."""
    design = np.asarray(design, dtype=float)
    coefficients = np.asarray(coefficients, dtype=float)
    logits = np.c_[design @ coefficients, np.zeros(len(design))]
    means = scipy.special.softmax(logits, axis=1)
    free_means = means[:, :-1]
    derivatives = np.asarray([
        np.diag(row) - np.outer(row, row) for row in free_means
    ])
    return means, derivatives


def _multinomial_gee_fit(
    counts,
    design,
    clusters,
    *,
    initial=None,
    max_iter=100,
    tolerance=1e-8,
    exchangeable=True,
):
    """Fit a multinomial marginal mean by Liang--Zeger GEE."""
    counts = np.asarray(counts, dtype=float)
    design = np.asarray(design, dtype=float)
    clusters = np.asarray(clusters)
    if counts.ndim != 2 or counts.shape[1] < 2:
        raise ValueError("counts must be samples by at least two paths")
    if (
        design.ndim != 2
        or design.shape[0] != len(counts)
        or len(clusters) != len(counts)
    ):
        raise ValueError("counts, design, and clusters must align")
    if np.any(counts < 0) or np.any(counts.sum(axis=1) <= 0):
        raise ValueError("each sample needs positive nonnegative counts")
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("GEE design is rank deficient")
    totals = counts.sum(axis=1)
    response = counts[:, :-1] / totals[:, None]
    n_paths = counts.shape[1]
    dimension = n_paths - 1
    n_coefficients = design.shape[1] * dimension
    cluster_rows = [
        np.flatnonzero(clusters == cluster)
        for cluster in np.unique(clusters)
    ]
    if len(cluster_rows) < 2:
        raise ValueError("GEE needs at least two clusters")
    if initial is None:
        fitted = _multinomial_fit(counts, design, max_iter=max_iter)
        coefficients = fitted["parameters"].reshape(
            design.shape[1], dimension
        )
    else:
        coefficients = np.asarray(initial, dtype=float).reshape(
            design.shape[1], dimension
        ).copy()
    alpha = 0.0
    scale = 1.0
    converged = False

    def components(current_coefficients, current_alpha, current_scale):
        means, mean_derivatives = _reference_multinomial_mean(
            design, current_coefficients
        )
        free_means = means[:, :-1]
        residuals = response - free_means
        observation_roots = []
        standardized = []
        derivative_rows = []
        for row in range(len(counts)):
            covariance = mean_derivatives[row] / totals[row]
            eigenvalues, eigenvectors = scipy.linalg.eigh(
                covariance, check_finite=False
            )
            floor = max(float(eigenvalues.max()) * 1e-10, 1e-12)
            eigenvalues = np.maximum(eigenvalues, floor)
            root = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T
            inverse_root = (
                eigenvectors * (1.0 / np.sqrt(eigenvalues))
            ) @ eigenvectors.T
            observation_roots.append(root)
            standardized.append(inverse_root @ residuals[row])
            derivative_rows.append(np.hstack([
                design[row, column] * mean_derivatives[row]
                for column in range(design.shape[1])
            ]))
        standardized = np.asarray(standardized)
        pearson_sum = float(np.sum(standardized * standardized))
        residual_df = max(
            len(counts) * dimension - n_coefficients, 1
        )
        estimated_scale = max(pearson_sum / residual_df, 1e-8)
        pair_sum = 0.0
        pair_count = 0
        for rows in cluster_rows:
            for first in range(len(rows)):
                for second in range(first):
                    pair_sum += float(
                        standardized[rows[first]]
                        @ standardized[rows[second]]
                    )
                    pair_count += dimension
        estimated_alpha = (
            pair_sum / (pair_count * estimated_scale)
            if pair_count
            else 0.0
        )
        maximum_cluster = max(len(rows) for rows in cluster_rows)
        lower = (
            -1.0 / (maximum_cluster - 1) + 1e-6
            if maximum_cluster > 1
            else -0.99
        )
        estimated_alpha = float(np.clip(estimated_alpha, lower, 0.95))
        use_alpha = estimated_alpha if exchangeable else 0.0
        use_scale = estimated_scale
        bread = np.zeros((n_coefficients, n_coefficients), dtype=float)
        score = np.zeros(n_coefficients, dtype=float)
        cluster_scores = []
        for rows in cluster_rows:
            size = len(rows)
            root = scipy.linalg.block_diag(*[
                observation_roots[row] for row in rows
            ])
            correlation = np.kron(
                np.full((size, size), use_alpha)
                + np.eye(size) * (1.0 - use_alpha),
                np.eye(dimension),
            )
            working = use_scale * root @ correlation @ root.T
            inverse = scipy.linalg.pinvh(
                working, rtol=1e-10, check_finite=False
            )
            derivative = np.vstack([
                derivative_rows[row] for row in rows
            ])
            residual = residuals[rows].ravel()
            cluster_score = derivative.T @ inverse @ residual
            bread += derivative.T @ inverse @ derivative
            score += cluster_score
            cluster_scores.append(cluster_score)
        return (
            means,
            estimated_alpha,
            estimated_scale,
            bread,
            score,
            cluster_scores,
        )

    iterations = 0
    for iterations in range(1, int(max_iter) + 1):
        (
            _,
            estimated_alpha,
            estimated_scale,
            bread,
            score,
            _,
        ) = components(coefficients, alpha, scale)
        step = scipy.linalg.pinvh(
            bread, rtol=1e-10, check_finite=False
        ) @ score
        coefficients += step.reshape(coefficients.shape)
        alpha = estimated_alpha
        scale = estimated_scale
        if np.linalg.norm(step, ord=np.inf) <= float(tolerance):
            converged = True
            break
    (
        means,
        alpha,
        scale,
        bread,
        score,
        cluster_scores,
    ) = components(coefficients, alpha, scale)
    bread_inverse = scipy.linalg.pinvh(
        bread, rtol=1e-10, check_finite=False
    )
    meat = sum(
        np.outer(cluster_score, cluster_score)
        for cluster_score in cluster_scores
    )
    robust_covariance = bread_inverse @ meat @ bread_inverse
    robust_covariance = 0.5 * (
        robust_covariance + robust_covariance.T
    )
    return {
        "coefficients": coefficients,
        "means": means,
        "working_correlation": alpha,
        "scale": scale,
        "model_covariance": bread_inverse,
        "robust_covariance": robust_covariance,
        "score_norm": float(np.linalg.norm(score, ord=np.inf)),
        "converged": bool(converged),
        "iterations": iterations,
        "clusters": len(cluster_rows),
    }


def multinomial_gee_test(
    counts,
    design,
    tested_columns,
    clusters,
    *,
    max_iter=100,
    tolerance=1e-8,
    exchangeable=True,
    initial=None,
):
    """Liang--Zeger multinomial GEE Wald test with sandwich covariance."""
    design = np.asarray(design, dtype=float)
    tested_columns = np.asarray(tested_columns, dtype=np.int64)
    fit = _multinomial_gee_fit(
        counts,
        design,
        clusters,
        max_iter=max_iter,
        tolerance=tolerance,
        exchangeable=exchangeable,
        initial=initial,
    )
    dimension = np.asarray(counts).shape[1] - 1
    tested = np.concatenate([
        np.arange(column * dimension, (column + 1) * dimension)
        for column in tested_columns
    ])
    values = fit["coefficients"].ravel()[tested]
    covariance = fit["robust_covariance"][np.ix_(tested, tested)]
    degrees_of_freedom = int(np.linalg.matrix_rank(covariance, tol=1e-10))
    if degrees_of_freedom:
        statistic = float(
            values @ scipy.linalg.pinvh(
                covariance, rtol=1e-10, check_finite=False
            ) @ values
        )
        p_value = float(
            scipy.stats.chi2.sf(statistic, degrees_of_freedom)
        )
    else:
        statistic = 0.0
        p_value = 1.0
    return {
        **fit,
        "statistic": statistic,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": p_value,
    }


def paired_logratio_test(
    counts,
    labels,
    clusters,
    *,
    pseudocount=0.5,
):
    """Test a paired compositional shift across independent subjects."""
    counts = np.asarray(counts, dtype=float)
    labels = np.asarray(labels)
    clusters = np.asarray(clusters)
    if counts.ndim != 2 or counts.shape[1] < 2:
        raise ValueError("counts must be samples by at least two categories")
    if len(counts) != len(labels) or len(counts) != len(clusters):
        raise ValueError("counts, labels, and clusters must align")
    levels = np.unique(labels)
    if len(levels) != 2:
        raise ValueError("paired log-ratio test requires two levels")
    differences = []
    for cluster in np.unique(clusters):
        rows = np.flatnonzero(clusters == cluster)
        first = rows[labels[rows] == levels[0]]
        second = rows[labels[rows] == levels[1]]
        if len(first) != 1 or len(second) != 1:
            continue
        if counts[first[0]].sum() <= 0 or counts[second[0]].sum() <= 0:
            continue
        logged = np.log(
            counts[[first[0], second[0]]] + float(pseudocount)
        )
        clr = logged - logged.mean(axis=1, keepdims=True)
        differences.append(clr[1] - clr[0])
    basis = scipy.linalg.helmert(counts.shape[1], full=False).T
    coordinates = (
        np.asarray(differences) @ basis
        if differences
        else np.empty((0, basis.shape[1]), dtype=float)
    )
    return paired_mean_test(coordinates)


def paired_mean_test(differences):
    """Test whether independent paired multivariate differences have mean zero."""
    differences = np.asarray(differences, dtype=float)
    if differences.ndim != 2 or differences.shape[1] < 1:
        raise ValueError("differences must be subjects by coordinates")
    n_subjects, dimension = differences.shape
    if n_subjects < 3 or n_subjects <= dimension:
        return {
            "statistic": 0.0,
            "degrees_of_freedom": dimension,
            "p_value": 1.0,
            "n_subjects": n_subjects,
            "converged": False,
        }
    mean = differences.mean(axis=0)
    if dimension == 1:
        standard_error = differences[:, 0].std(ddof=1) / np.sqrt(n_subjects)
        statistic = (
            float(mean[0] / standard_error)
            if standard_error > 0
            else 0.0
        )
        p_value = float(
            2.0 * scipy.stats.t.sf(abs(statistic), n_subjects - 1)
        )
        reported = statistic * statistic
    else:
        covariance = np.cov(differences, rowvar=False, ddof=1)
        if np.linalg.matrix_rank(covariance, tol=1e-10) < dimension:
            return {
                "statistic": 0.0,
                "degrees_of_freedom": dimension,
                "p_value": 1.0,
                "n_subjects": n_subjects,
                "converged": False,
            }
        t_squared = float(
            n_subjects
            * mean
            @ scipy.linalg.solve(
                covariance, mean, assume_a="sym", check_finite=False
            )
        )
        statistic = (
            (n_subjects - dimension)
            * t_squared
            / (dimension * (n_subjects - 1))
        )
        p_value = float(
            scipy.stats.f.sf(
                statistic, dimension, n_subjects - dimension
            )
        )
        reported = t_squared
    return {
        "statistic": reported,
        "degrees_of_freedom": dimension,
        "p_value": p_value,
        "n_subjects": n_subjects,
        "converged": True,
    }


def fit_variance_prior(variances, residual_df, *, winsor_tail=0.05):
    """Fit a scaled inverse-chi-squared prior to ``M`` sample variances.

    ``variances`` and ``residual_df`` are length-``M`` vectors. The method
    matches the first two log-variance moments after removing each test's
    chi-squared sampling moments. Symmetric Winsorization limits the influence
    of genes with unusual biological variance.
    """
    variances = np.asarray(variances, dtype=float)
    residual_df = np.broadcast_to(
        np.asarray(residual_df, dtype=float), variances.shape
    )
    valid = (
        np.isfinite(variances)
        & (variances > 0)
        & np.isfinite(residual_df)
        & (residual_df > 0)
    )
    if valid.sum() < 3:
        raise ValueError("at least three finite positive variances are required")
    variances = variances[valid]
    residual_df = residual_df[valid]
    winsor_tail = float(winsor_tail)
    if not 0 <= winsor_tail < 0.5:
        raise ValueError("winsor_tail must be in [0, 0.5)")
    adjusted = (
        np.log(variances)
        - scipy.special.digamma(residual_df / 2)
        + np.log(residual_df / 2)
    )
    if winsor_tail:
        lower, upper = np.quantile(
            adjusted, [winsor_tail, 1 - winsor_tail]
        )
        adjusted = np.clip(adjusted, lower, upper)
    target = max(
        float(
            adjusted.var(ddof=1)
            - np.mean(scipy.special.polygamma(1, residual_df / 2))
        ),
        1e-8,
    )
    prior_df = float(
        scipy.optimize.brentq(
            lambda value: scipy.special.polygamma(1, value / 2) - target,
            1e-3,
            1e12,
        )
    )
    prior_variance = float(
        np.exp(
            adjusted.mean()
            + scipy.special.digamma(prior_df / 2)
            - np.log(prior_df / 2)
        )
    )
    return prior_df, prior_variance


def moderated_t_pvalues(
    t_statistics,
    variances,
    residual_df,
    prior_df,
    prior_variance,
):
    """Return two-sided moderated-t p-values for ``M`` scalar tests."""
    t_statistics = np.asarray(t_statistics, dtype=float)
    variances = np.broadcast_to(
        np.asarray(variances, dtype=float), t_statistics.shape
    )
    residual_df = np.broadcast_to(
        np.asarray(residual_df, dtype=float), t_statistics.shape
    )
    prior_df = float(prior_df)
    prior_variance = float(prior_variance)
    if prior_df <= 0 or prior_variance <= 0:
        raise ValueError("prior_df and prior_variance must be positive")
    posterior = (
        residual_df * variances + prior_df * prior_variance
    ) / (residual_df + prior_df)
    moderated = np.abs(t_statistics) * np.sqrt(variances / posterior)
    return 2 * scipy.stats.t.sf(moderated, residual_df + prior_df)


def _multinomial_glmm_cluster_mode(
    counts,
    fixed_logits,
    variance,
    *,
    initial=None,
    max_iter=50,
    tolerance=1e-10,
):
    """Find one cluster's Gaussian random-intercept posterior mode."""
    counts = np.asarray(counts, dtype=float)
    fixed_logits = np.asarray(fixed_logits, dtype=float)
    dimension = counts.shape[1] - 1
    totals = counts.sum(axis=1)
    mode = (
        np.zeros(dimension, dtype=float)
        if initial is None
        else np.asarray(initial, dtype=float).copy()
    )
    precision = 1.0 / float(variance)
    hessian = np.eye(dimension) * precision

    def log_posterior(value):
        logits = np.c_[
            fixed_logits + value[None, :],
            np.zeros(len(counts)),
        ]
        log_means = logits - scipy.special.logsumexp(
            logits, axis=1, keepdims=True
        )
        return (
            float(np.sum(counts * log_means))
            - 0.5 * precision * float(value @ value)
        )

    for _ in range(int(max_iter)):
        logits = np.c_[
            fixed_logits + mode[None, :],
            np.zeros(len(counts)),
        ]
        means = scipy.special.softmax(logits, axis=1)
        free = means[:, :-1]
        gradient = (
            counts[:, :-1] - totals[:, None] * free
        ).sum(axis=0) - precision * mode
        hessian = np.eye(dimension) * precision
        for total, row in zip(totals, free):
            hessian += total * (
                np.diag(row) - np.outer(row, row)
            )
        step = scipy.linalg.solve(
            hessian, gradient, assume_a="pos", check_finite=False
        )
        current = log_posterior(mode)
        directional = float(gradient @ step)
        step_size = 1.0
        while step_size > 1e-8:
            candidate = mode + step_size * step
            if log_posterior(candidate) >= (
                current + 1e-4 * step_size * directional
            ):
                break
            step_size *= 0.5
        mode += step_size * step
        if (
            np.linalg.norm(step_size * step, ord=np.inf)
            <= float(tolerance)
        ):
            break
    logits = np.c_[fixed_logits + mode[None, :], np.zeros(len(counts))]
    means = scipy.special.softmax(logits, axis=1)[:, :-1]
    hessian = np.eye(dimension) * precision
    for total, row in zip(totals, means):
        hessian += total * (np.diag(row) - np.outer(row, row))
    return mode, hessian


_GLMM_JAX_FUNCTIONS = {}


def _multinomial_glmm_jax_functions(
    n_samples, n_paths, n_design_columns, n_clusters
):
    """Build and cache shape-specialized JAX Laplace functions."""
    key = (n_samples, n_paths, n_design_columns, n_clusters)
    if key in _GLMM_JAX_FUNCTIONS:
        return _GLMM_JAX_FUNCTIONS[key]
    import jax
    import jax.numpy as jnp
    import jax.scipy as jsp

    jax.config.update("jax_enable_x64", True)
    dimension = n_paths - 1
    coefficient_count = n_design_columns * dimension
    identity = jnp.eye(dimension, dtype=jnp.float64)

    def posterior_modes(
        parameters, counts_jax, design_jax, cluster_jax, membership
    ):
        totals = counts_jax.sum(axis=1)
        coefficients = parameters[:coefficient_count].reshape(
            n_design_columns, dimension
        )
        variance = jnp.exp(2.0 * parameters[-1])
        fixed_logits = design_jax @ coefficients
        modes = jnp.zeros((n_clusters, dimension), dtype=jnp.float64)

        def update(_, current_modes):
            logits = jnp.concatenate(
                (
                    fixed_logits + current_modes[cluster_jax],
                    jnp.zeros((n_samples, 1), dtype=jnp.float64),
                ),
                axis=1,
            )
            means = jax.nn.softmax(logits, axis=1)[:, :-1]
            residual = counts_jax[:, :-1] - totals[:, None] * means
            gradient = membership.T @ residual - current_modes / variance
            observation_hessian = jax.vmap(
                lambda total, row: total
                * (jnp.diag(row) - jnp.outer(row, row))
            )(totals, means)
            hessian = (
                jnp.einsum("im,ijk->mjk", membership, observation_hessian)
                + identity[None, :, :] / variance
            )
            step = jnp.linalg.solve(hessian, gradient[..., None])[..., 0]
            maximum = jnp.max(jnp.abs(step), axis=1, keepdims=True)
            damping = jnp.minimum(1.0, 2.0 / (maximum + 1e-12))
            return current_modes + damping * step

        modes = jax.lax.fori_loop(0, 30, update, modes)
        return fixed_logits, modes, variance

    def objective(
        parameters, counts_jax, design_jax, cluster_jax, membership
    ):
        totals = counts_jax.sum(axis=1)
        fixed_logits, modes, variance = posterior_modes(
            parameters, counts_jax, design_jax, cluster_jax, membership
        )
        logits = jnp.concatenate(
            (
                fixed_logits + modes[cluster_jax],
                jnp.zeros((n_samples, 1), dtype=jnp.float64),
            ),
            axis=1,
        )
        log_means = logits - jsp.special.logsumexp(
            logits, axis=1, keepdims=True
        )
        means = jnp.exp(log_means)[:, :-1]
        observation_hessian = jax.vmap(
            lambda total, row: total
            * (jnp.diag(row) - jnp.outer(row, row))
        )(totals, means)
        hessian = (
            jnp.einsum("im,ijk->mjk", membership, observation_hessian)
            + identity[None, :, :] / variance
        )
        log_determinant = jnp.linalg.slogdet(hessian)[1]
        value = (
            jnp.sum(counts_jax * log_means)
            - 0.5 * jnp.sum(modes * modes) / variance
            - n_clusters * dimension * parameters[-1]
            - 0.5 * jnp.sum(log_determinant)
        )
        return -value

    functions = (
        jax.jit(jax.value_and_grad(objective)),
        jax.jit(posterior_modes),
        jnp,
        jax,
    )
    _GLMM_JAX_FUNCTIONS[key] = functions
    return functions


def _multinomial_glmm_jax_fit(
    counts,
    design,
    cluster_index,
    n_clusters,
    initial,
    max_iter,
):
    """Autodifferentiated Laplace fit used when JAX is available."""
    value_and_gradient, posterior_modes, jnp, jax = (
        _multinomial_glmm_jax_functions(
            len(counts),
            counts.shape[1],
            design.shape[1],
            int(n_clusters),
        )
    )
    counts_jax = jnp.asarray(counts, dtype=jnp.float64)
    design_jax = jnp.asarray(design, dtype=jnp.float64)
    cluster_jax = jnp.asarray(cluster_index, dtype=jnp.int32)
    membership = jax.nn.one_hot(
        cluster_jax, int(n_clusters), dtype=jnp.float64
    )
    dimension = counts.shape[1] - 1
    coefficient_count = design.shape[1] * dimension

    def scipy_objective(parameters):
        value, gradient = value_and_gradient(
            jnp.asarray(parameters),
            counts_jax,
            design_jax,
            cluster_jax,
            membership,
        )
        return float(value), np.asarray(gradient, dtype=float)

    bounds = [(-20.0, 20.0)] * coefficient_count + [(-8.0, 4.0)]
    result = scipy.optimize.minimize(
        scipy_objective,
        np.asarray(initial, dtype=float),
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": int(max_iter), "ftol": 1e-9, "gtol": 1e-5},
    )
    final_value, final_gradient = scipy_objective(result.x)
    fixed_logits, modes, variance = posterior_modes(
        jnp.asarray(result.x),
        counts_jax,
        design_jax,
        cluster_jax,
        membership,
    )
    fixed_logits = np.asarray(fixed_logits)
    modes = np.asarray(modes)
    variance = float(variance)
    mode_score = 0.0
    for cluster in range(int(n_clusters)):
        rows = np.flatnonzero(cluster_index == cluster)
        logits = np.c_[
            fixed_logits[rows] + modes[cluster][None, :],
            np.zeros(len(rows)),
        ]
        means = scipy.special.softmax(logits, axis=1)[:, :-1]
        score = (
            counts[rows, :-1]
            - counts[rows].sum(axis=1)[:, None] * means
        ).sum(axis=0) - modes[cluster] / variance
        mode_score = max(
            mode_score, float(np.linalg.norm(score, ord=np.inf))
        )
    gradient_norm = float(np.linalg.norm(final_gradient, ord=np.inf))
    return {
        "objective": final_value,
        "parameters": np.asarray(result.x),
        "coefficients": np.asarray(result.x[:coefficient_count]).reshape(
            design.shape[1], dimension
        ),
        "random_effect_sd": float(np.sqrt(variance)),
        "random_effect_modes": modes,
        "mode_score_norm": mode_score,
        "gradient_norm": gradient_norm,
        "converged": bool(
            result.success
            or (
                np.isfinite(final_value)
                and gradient_norm <= 1e-3
                and mode_score <= 1e-6
            )
        ),
        "iterations": int(result.nit),
        "message": str(result.message),
    }


def _multinomial_glmm_laplace_fit(
    counts,
    design,
    clusters,
    *,
    initial=None,
    max_iter=200,
):
    """Fit a multinomial-logit random-intercept GLMM by Laplace ML."""
    counts = np.asarray(counts, dtype=float)
    design = np.asarray(design, dtype=float)
    clusters = np.asarray(clusters)
    if counts.ndim != 2 or counts.shape[1] < 2:
        raise ValueError("counts must be samples by at least two paths")
    if (
        design.ndim != 2
        or design.shape[0] != len(counts)
        or len(clusters) != len(counts)
    ):
        raise ValueError("counts, design, and clusters must align")
    if np.any(counts < 0) or np.any(counts.sum(axis=1) <= 0):
        raise ValueError("each sample needs positive nonnegative counts")
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("GLMM design is rank deficient")
    unique_clusters, cluster_index = np.unique(
        clusters, return_inverse=True
    )
    cluster_rows = [
        np.flatnonzero(cluster_index == cluster)
        for cluster in range(len(unique_clusters))
    ]
    dimension = counts.shape[1] - 1
    coefficient_count = design.shape[1] * dimension
    if initial is None:
        fixed = _multinomial_fit(counts, design, max_iter=max_iter)
        parameters = np.r_[fixed["parameters"], np.log(0.5)]
    else:
        parameters = np.asarray(initial, dtype=float).copy()
    if len(parameters) != coefficient_count + 1:
        raise ValueError("initial GLMM parameters have the wrong length")
    try:
        return _multinomial_glmm_jax_fit(
            counts,
            design,
            cluster_index,
            len(unique_clusters),
            parameters,
            max_iter,
        )
    except ImportError:
        pass
    cached_modes = np.zeros((len(cluster_rows), dimension), dtype=float)

    def objective(current):
        coefficients = current[:coefficient_count].reshape(
            design.shape[1], dimension
        )
        log_sd = float(current[-1])
        variance = float(np.exp(2.0 * log_sd))
        fixed_logits = design @ coefficients
        value = 0.0
        modes = np.zeros_like(cached_modes)
        maximum_mode_score = 0.0
        for cluster, rows in enumerate(cluster_rows):
            mode, hessian = _multinomial_glmm_cluster_mode(
                counts[rows],
                fixed_logits[rows],
                variance,
                initial=cached_modes[cluster],
            )
            modes[cluster] = mode
            logits = np.c_[
                fixed_logits[rows] + mode[None, :],
                np.zeros(len(rows)),
            ]
            log_means = logits - scipy.special.logsumexp(
                logits, axis=1, keepdims=True
            )
            log_likelihood = float(np.sum(counts[rows] * log_means))
            sign, log_determinant = np.linalg.slogdet(hessian)
            if sign <= 0:
                return np.inf
            value += (
                log_likelihood
                - 0.5 * np.dot(mode, mode) / variance
                - dimension * log_sd
                - 0.5 * log_determinant
            )
            totals = counts[rows].sum(axis=1)
            means = np.exp(log_means)[:, :-1]
            score = (
                counts[rows, :-1] - totals[:, None] * means
            ).sum(axis=0) - mode / variance
            maximum_mode_score = max(
                maximum_mode_score,
                float(np.linalg.norm(score, ord=np.inf)),
            )
        cached_modes[:] = modes
        objective.last_modes = modes
        objective.last_mode_score = maximum_mode_score
        return -value

    bounds = [(-20.0, 20.0)] * coefficient_count + [(-8.0, 4.0)]
    result = scipy.optimize.minimize(
        objective,
        parameters,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": int(max_iter), "ftol": 1e-9, "gtol": 1e-5},
    )
    final_objective = objective(result.x)
    gradient_norm = (
        float(np.linalg.norm(result.jac, ord=np.inf))
        if result.jac is not None
        else np.inf
    )
    return {
        "objective": final_objective,
        "parameters": np.asarray(result.x),
        "coefficients": np.asarray(result.x[:coefficient_count]).reshape(
            design.shape[1], dimension
        ),
        "random_effect_sd": float(np.exp(result.x[-1])),
        "random_effect_modes": np.asarray(objective.last_modes),
        "mode_score_norm": float(objective.last_mode_score),
        "gradient_norm": gradient_norm,
        "converged": bool(
            result.success
            or (
                np.isfinite(final_objective)
                and gradient_norm <= 1e-3
                and objective.last_mode_score <= 1e-6
            )
        ),
        "iterations": int(result.nit),
        "message": str(result.message),
    }


def multinomial_glmm_test(
    counts,
    null_design,
    alternative_design,
    clusters,
    *,
    max_iter=200,
    fitted_null=None,
):
    """Laplace likelihood-ratio test for a multinomial mouse-intercept GLMM."""
    counts = np.asarray(counts, dtype=float)
    null_design = np.asarray(null_design, dtype=float)
    alternative_design = np.asarray(alternative_design, dtype=float)
    if fitted_null is None:
        null = _multinomial_glmm_laplace_fit(
            counts, null_design, clusters, max_iter=max_iter
        )
    else:
        null = {
            "objective": float(fitted_null["null_objective"]),
            "coefficients": np.asarray(
                fitted_null["null_coefficients"], dtype=float
            ),
            "parameters": np.r_[
                np.asarray(
                    fitted_null["null_coefficients"], dtype=float
                ).ravel(),
                np.log(float(fitted_null["null_random_effect_sd"])),
            ],
            "random_effect_sd": float(
                fitted_null["null_random_effect_sd"]
            ),
            "converged": bool(fitted_null["null_converged"]),
            "gradient_norm": float(
                fitted_null.get("null_gradient_norm", np.nan)
            ),
            "mode_score_norm": float(
                fitted_null.get("null_mode_score_norm", np.nan)
            ),
            "iterations": 0,
        }
    dimension = counts.shape[1] - 1
    initial_coefficients = np.zeros(
        (alternative_design.shape[1], dimension), dtype=float
    )
    initial_coefficients[: null_design.shape[1]] = null["coefficients"]
    initial = np.r_[initial_coefficients.ravel(), null["parameters"][-1]]
    alternative = _multinomial_glmm_laplace_fit(
        counts,
        alternative_design,
        clusters,
        initial=initial,
        max_iter=max_iter,
    )
    statistic = max(
        2.0 * (null["objective"] - alternative["objective"]), 0.0
    )
    degrees_of_freedom = (
        alternative_design.shape[1] - null_design.shape[1]
    ) * dimension
    return {
        "statistic": statistic,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": float(
            scipy.stats.chi2.sf(statistic, degrees_of_freedom)
        ),
        "null_converged": null["converged"],
        "alternative_converged": alternative["converged"],
        "null_objective": null["objective"],
        "alternative_objective": alternative["objective"],
        "null_random_effect_sd": null["random_effect_sd"],
        "alternative_random_effect_sd": alternative["random_effect_sd"],
        "null_coefficients": null["coefficients"],
        "alternative_coefficients": alternative["coefficients"],
        "null_gradient_norm": null["gradient_norm"],
        "alternative_gradient_norm": alternative["gradient_norm"],
        "null_mode_score_norm": null["mode_score_norm"],
        "alternative_mode_score_norm": alternative["mode_score_norm"],
        "null_iterations": null["iterations"],
        "alternative_iterations": alternative["iterations"],
    }


def effective_multinomial_size(proportions, covariance, maximum=None):
    """Match ILR covariance to a multinomial and return its effective size."""
    proportions = np.asarray(proportions, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if (
        proportions.ndim != 1
        or len(proportions) < 2
        or np.any(proportions <= 0)
        or not np.isclose(proportions.sum(), 1.0)
    ):
        raise ValueError("proportions must be a positive probability vector")
    basis = helmert_basis(len(proportions))
    unit_covariance = basis.T @ np.diag(1.0 / proportions) @ basis
    if covariance.shape != unit_covariance.shape:
        raise ValueError("covariance has the wrong shape")
    inverse_size = float(
        np.sum(unit_covariance * covariance)
        / np.sum(unit_covariance * unit_covariance)
    )
    if not np.isfinite(inverse_size) or inverse_size <= 0:
        raise ValueError("covariance does not imply a positive effective size")
    size = 1.0 / inverse_size
    if maximum is not None:
        size = min(size, float(maximum))
    return size


def integerize_compositional_counts(counts):
    """Round fractional compositional counts by the largest-remainder rule."""
    counts = np.asarray(counts, dtype=float)
    if (
        counts.ndim != 2
        or counts.shape[1] < 2
        or np.any(counts < 0)
        or np.any(counts.sum(axis=1) <= 0)
    ):
        raise ValueError(
            "counts must be samples by paths with positive row totals"
        )
    integer_totals = np.maximum(
        np.rint(counts.sum(axis=1)).astype(np.int64), 1
    )
    proportions = counts / counts.sum(axis=1, keepdims=True)
    targets = integer_totals[:, None] * proportions
    result = np.floor(targets).astype(np.int64)
    for row in range(len(result)):
        remainder = int(integer_totals[row] - result[row].sum())
        if remainder:
            order = np.argsort(
                -(targets[row] - result[row]), kind="stable"
            )
            result[row, order[:remainder]] += 1
    return result


def project_composition(proportions, covariance, selected):
    """Project an ILR composition and covariance onto selected categories."""
    proportions = np.asarray(proportions, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    selected = np.asarray(selected, dtype=np.int64)
    if (
        proportions.ndim != 1
        or len(proportions) < 2
        or covariance.shape != (len(proportions) - 1,) * 2
    ):
        raise ValueError("composition and ILR covariance shapes disagree")
    if (
        selected.ndim != 1
        or len(selected) < 2
        or len(np.unique(selected)) != len(selected)
        or np.any(selected < 0)
        or np.any(selected >= len(proportions))
    ):
        raise ValueError("selected categories must be unique valid indices")
    old_basis = helmert_basis(len(proportions))
    new_basis = helmert_basis(len(selected))
    embedding = np.eye(len(proportions))[:, selected]
    transform = new_basis.T @ embedding.T @ old_basis
    projected_covariance = transform @ covariance @ transform.T
    projected_proportions = proportions[selected].copy()
    projected_proportions /= projected_proportions.sum()
    return projected_proportions, projected_covariance


def permutation_rank_p_value(
    observed_statistic,
    null_statistics,
    *,
    relative_tolerance=1e-6,
):
    """Return a finite-permutation p-value with stable numerical ties."""
    observed_statistic = float(observed_statistic)
    null_statistics = np.asarray(null_statistics, dtype=float)
    if null_statistics.ndim != 1 or np.any(np.isnan(null_statistics)):
        raise ValueError("null statistics must be a one-dimensional array")
    tolerance = float(relative_tolerance) * max(
        1.0, abs(observed_statistic)
    )
    return float(
        (1 + np.sum(null_statistics >= observed_statistic - tolerance))
        / (len(null_statistics) + 1)
    )


def _dirichlet_multinomial_fit(
    counts,
    design,
    initial=None,
    max_iter=500,
    fixed_concentration=None,
):
    counts = np.asarray(counts, dtype=float)
    design = np.asarray(design, dtype=float)
    if counts.ndim != 2 or counts.shape[1] < 2:
        raise ValueError("counts must be samples by at least two paths")
    if design.ndim != 2 or design.shape[0] != counts.shape[0]:
        raise ValueError("design and counts have different sample counts")
    if np.any(counts < 0) or np.any(counts.sum(axis=1) <= 0):
        raise ValueError("each sample needs positive nonnegative counts")
    n_paths = counts.shape[1]
    n_coefficients = design.shape[1] * (n_paths - 1)
    fit_concentration = fixed_concentration is None
    if initial is None:
        pooled = counts.sum(axis=0) + 0.5
        pooled /= pooled.sum()
        coefficients = np.zeros((design.shape[1], n_paths - 1))
        if np.allclose(design[:, 0], 1.0):
            coefficients[0] = np.log(pooled[:-1] / pooled[-1])
        initial = coefficients.ravel()
        if fit_concentration:
            initial = np.r_[initial, np.log(10.0)]
    else:
        initial = np.asarray(initial, dtype=float)
    expected_parameters = n_coefficients + int(fit_concentration)
    if len(initial) != expected_parameters:
        raise ValueError("initial value has the wrong length")
    totals = counts.sum(axis=1)

    def objective(parameters):
        coefficients = parameters[:n_coefficients].reshape(
            design.shape[1], n_paths - 1
        )
        logits = np.c_[design @ coefficients, np.zeros(len(counts))]
        means = scipy.special.softmax(logits, axis=1)
        concentration = (
            np.exp(parameters[-1])
            if fit_concentration
            else float(fixed_concentration)
        )
        alpha = concentration * means
        log_likelihood = np.sum(
            scipy.special.gammaln(concentration)
            - scipy.special.gammaln(concentration + totals)
            + np.sum(
                scipy.special.gammaln(alpha + counts)
                - scipy.special.gammaln(alpha),
                axis=1,
            )
        )
        alpha_score = concentration * (
            scipy.special.digamma(alpha + counts)
            - scipy.special.digamma(alpha)
        )
        centered_score = alpha_score - np.sum(
            means * alpha_score, axis=1, keepdims=True
        )
        logit_score = means * centered_score
        coefficient_score = design.T @ logit_score[:, :-1]
        gradient = coefficient_score.ravel()
        if fit_concentration:
            concentration_score = concentration * np.sum(
                scipy.special.digamma(concentration)
                - scipy.special.digamma(concentration + totals)
                + np.sum(
                    means
                    * (
                        scipy.special.digamma(alpha + counts)
                        - scipy.special.digamma(alpha)
                    ),
                    axis=1,
                )
            )
            gradient = np.r_[gradient, concentration_score]
        return -float(log_likelihood), -gradient

    bounds = [(-20.0, 20.0)] * n_coefficients
    if fit_concentration:
        bounds.append((-10.0, np.log(1e6)))
    result = scipy.optimize.minimize(
        objective,
        initial,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": int(max_iter), "ftol": 1e-10},
    )
    gradient_norm = float(np.linalg.norm(result.jac, ord=np.inf))
    return {
        "objective": float(result.fun),
        "parameters": np.asarray(result.x),
        "concentration": (
            float(np.exp(result.x[-1]))
            if fit_concentration
            else float(fixed_concentration)
        ),
        "converged": bool(result.success or gradient_norm <= 1e-4),
        "iterations": int(result.nit),
        "message": str(result.message),
        "gradient_norm": gradient_norm,
    }


def _multinomial_fit(counts, design, initial=None, max_iter=500):
    counts = np.asarray(counts, dtype=float)
    design = np.asarray(design, dtype=float)
    n_paths = counts.shape[1]
    n_coefficients = design.shape[1] * (n_paths - 1)
    if initial is None:
        pooled = counts.sum(axis=0) + 0.5
        pooled /= pooled.sum()
        coefficients = np.zeros((design.shape[1], n_paths - 1))
        if np.allclose(design[:, 0], 1.0):
            coefficients[0] = np.log(pooled[:-1] / pooled[-1])
        initial = coefficients.ravel()
    initial = np.asarray(initial, dtype=float)
    totals = counts.sum(axis=1)

    def objective(parameters):
        coefficients = parameters.reshape(design.shape[1], n_paths - 1)
        logits = np.c_[design @ coefficients, np.zeros(len(counts))]
        log_means = logits - scipy.special.logsumexp(
            logits, axis=1, keepdims=True
        )
        means = np.exp(log_means)
        log_likelihood = float(np.sum(counts * log_means))
        logit_score = counts - totals[:, None] * means
        gradient = design.T @ logit_score[:, :-1]
        return -log_likelihood, -gradient.ravel()

    result = scipy.optimize.minimize(
        objective,
        initial,
        method="L-BFGS-B",
        jac=True,
        bounds=[(-20.0, 20.0)] * n_coefficients,
        options={"maxiter": int(max_iter), "ftol": 1e-10},
    )
    gradient_norm = float(np.linalg.norm(result.jac, ord=np.inf))
    return {
        "objective": float(result.fun),
        "parameters": np.asarray(result.x),
        "converged": bool(result.success or gradient_norm <= 1e-4),
        "iterations": int(result.nit),
        "message": str(result.message),
        "gradient_norm": gradient_norm,
    }


def multinomial_test(
    counts,
    null_design,
    alternative_design,
    *,
    max_iter=500,
    fitted_null=None,
):
    """Likelihood-ratio test for multinomial compositional regression."""
    counts = np.asarray(counts, dtype=float)
    null_design = np.asarray(null_design, dtype=float)
    alternative_design = np.asarray(alternative_design, dtype=float)
    if (
        null_design.shape[0] != len(counts)
        or alternative_design.shape[0] != len(counts)
    ):
        raise ValueError("design and counts have different sample counts")
    if np.linalg.matrix_rank(null_design) != null_design.shape[1]:
        raise ValueError("null design is rank deficient")
    if np.linalg.matrix_rank(alternative_design) != alternative_design.shape[1]:
        raise ValueError("alternative design is rank deficient")
    if (
        alternative_design.shape[1] < null_design.shape[1]
        or not np.allclose(
            alternative_design[:, : null_design.shape[1]],
            null_design,
        )
    ):
        raise ValueError("alternative design must begin with the null design")
    if fitted_null is None:
        null = _multinomial_fit(
            counts,
            null_design,
            max_iter=max_iter,
        )
        if not null["converged"]:
            null = _multinomial_fit(
                counts,
                null_design,
                initial=null["parameters"],
                max_iter=4 * int(max_iter),
            )
    else:
        null_coefficients = np.asarray(
            fitted_null["null_coefficients"], dtype=float
        )
        if null_coefficients.shape != (
            null_design.shape[1], counts.shape[1] - 1
        ):
            raise ValueError("fitted null coefficients have the wrong shape")
        null = {
            "objective": float(fitted_null["null_objective"]),
            "parameters": null_coefficients.ravel(),
            "converged": True,
            "iterations": 0,
            "message": "reused fitted null",
        }
    n_paths = counts.shape[1]
    null_coefficients = null["parameters"].reshape(
        null_design.shape[1], n_paths - 1
    )
    initial = np.zeros((alternative_design.shape[1], n_paths - 1))
    initial[: len(null_coefficients)] = null_coefficients
    alternative = _multinomial_fit(
        counts,
        alternative_design,
        initial=initial.ravel(),
        max_iter=max_iter,
    )
    if not alternative["converged"]:
        alternative = _multinomial_fit(
            counts,
            alternative_design,
            initial=alternative["parameters"],
            max_iter=4 * int(max_iter),
        )
    statistic = max(
        2.0 * (null["objective"] - alternative["objective"]),
        0.0,
    )
    degrees_of_freedom = (
        alternative_design.shape[1] - null_design.shape[1]
    ) * (n_paths - 1)
    return {
        "model": "multinomial",
        "statistic": statistic,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": float(scipy.stats.chi2.sf(statistic, degrees_of_freedom)),
        "null_converged": null["converged"],
        "alternative_converged": alternative["converged"],
        "null_iterations": null["iterations"],
        "alternative_iterations": alternative["iterations"],
        "null_message": null["message"],
        "alternative_message": alternative["message"],
        "null_objective": null["objective"],
        "null_coefficients": null["parameters"].reshape(
            null_design.shape[1], n_paths - 1
        ),
        "alternative_coefficients": alternative["parameters"].reshape(
            alternative_design.shape[1], n_paths - 1
        ),
    }


def dirichlet_multinomial_test(
    counts,
    null_design,
    alternative_design,
    *,
    max_iter=500,
    fix_null_concentration=True,
    fitted_null=None,
):
    """Likelihood-ratio test of compositional regression coefficients."""
    counts = np.asarray(counts, dtype=float)
    null_design = np.asarray(null_design, dtype=float)
    alternative_design = np.asarray(alternative_design, dtype=float)
    if (
        null_design.shape[0] != len(counts)
        or alternative_design.shape[0] != len(counts)
    ):
        raise ValueError("design and counts have different sample counts")
    if np.linalg.matrix_rank(null_design) != null_design.shape[1]:
        raise ValueError("null design is rank deficient")
    if np.linalg.matrix_rank(alternative_design) != alternative_design.shape[1]:
        raise ValueError("alternative design is rank deficient")
    if (
        alternative_design.shape[1] < null_design.shape[1]
        or not np.allclose(
            alternative_design[:, : null_design.shape[1]],
            null_design,
        )
    ):
        raise ValueError("alternative design must begin with the null design")
    if fitted_null is None:
        null = _dirichlet_multinomial_fit(
            counts,
            null_design,
            max_iter=max_iter,
        )
        if not null["converged"]:
            null = _dirichlet_multinomial_fit(
                counts,
                null_design,
                initial=null["parameters"],
                max_iter=4 * int(max_iter),
            )
    elif fitted_null["model"] == "multinomial":
        result = multinomial_test(
            counts,
            null_design,
            alternative_design,
            max_iter=max_iter,
            fitted_null=fitted_null,
        )
        result.update({
            "null_concentration": np.inf,
            "alternative_concentration": np.inf,
        })
        return result
    else:
        null_coefficients = np.asarray(
            fitted_null["null_coefficients"], dtype=float
        )
        concentration = float(fitted_null["null_concentration"])
        if null_coefficients.shape != (
            null_design.shape[1], counts.shape[1] - 1
        ):
            raise ValueError("fitted null coefficients have the wrong shape")
        null = {
            "objective": float(fitted_null["null_objective"]),
            "parameters": np.r_[
                null_coefficients.ravel(),
                np.log(min(concentration, 1e6)),
            ],
            "concentration": concentration,
            "converged": True,
            "iterations": 0,
            "message": "reused fitted null",
        }
    if (
        null["concentration"] >= 0.99e6
        and fix_null_concentration
    ):
        result = multinomial_test(
            counts,
            null_design,
            alternative_design,
            max_iter=max_iter,
        )
        result.update({
            "null_concentration": np.inf,
            "alternative_concentration": np.inf,
        })
        return result
    if null["concentration"] >= 0.99e6:
        multinomial_null = _multinomial_fit(
            counts,
            null_design,
            initial=null["parameters"][:-1],
            max_iter=max_iter,
        )
        if not multinomial_null["converged"]:
            multinomial_null = _multinomial_fit(
                counts,
                null_design,
                initial=multinomial_null["parameters"],
                max_iter=4 * int(max_iter),
            )
        null = {
            **multinomial_null,
            "parameters": np.r_[
                multinomial_null["parameters"], np.log(1e6)
            ],
            "concentration": np.inf,
        }
    n_paths = counts.shape[1]
    null_coefficients = null["parameters"][:-1].reshape(
        null_design.shape[1], n_paths - 1
    )
    initial_coefficients = np.zeros(
        (alternative_design.shape[1], n_paths - 1)
    )
    initial_coefficients[: len(null_coefficients)] = null_coefficients
    initial = initial_coefficients.ravel()
    if not fix_null_concentration:
        initial = np.r_[initial, null["parameters"][-1]]
    alternative = _dirichlet_multinomial_fit(
        counts,
        alternative_design,
        initial=initial,
        max_iter=max_iter,
        fixed_concentration=(
            null["concentration"] if fix_null_concentration else None
        ),
    )
    if not alternative["converged"]:
        alternative = _dirichlet_multinomial_fit(
            counts,
            alternative_design,
            initial=alternative["parameters"],
            max_iter=4 * int(max_iter),
            fixed_concentration=(
                null["concentration"]
                if fix_null_concentration
                else None
            ),
        )
    if not fix_null_concentration:
        multinomial_alternative = _multinomial_fit(
            counts,
            alternative_design,
            initial=alternative["parameters"][:
                alternative_design.shape[1] * (n_paths - 1)
            ],
            max_iter=max_iter,
        )
        if not multinomial_alternative["converged"]:
            multinomial_alternative = _multinomial_fit(
                counts,
                alternative_design,
                initial=multinomial_alternative["parameters"],
                max_iter=4 * int(max_iter),
            )
        if (
            multinomial_alternative["converged"]
            and (
                not alternative["converged"]
                or multinomial_alternative["objective"]
                <= alternative["objective"]
            )
        ):
            alternative = {
                **multinomial_alternative,
                "concentration": np.inf,
            }
    statistic = max(
        2.0 * (null["objective"] - alternative["objective"]),
        0.0,
    )
    degrees_of_freedom = (
        alternative_design.shape[1] - null_design.shape[1]
    ) * (n_paths - 1)
    return {
        "model": "dirichlet_multinomial",
        "statistic": statistic,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": float(scipy.stats.chi2.sf(statistic, degrees_of_freedom)),
        "null_concentration": null["concentration"],
        "alternative_concentration": alternative["concentration"],
        "null_converged": null["converged"],
        "alternative_converged": alternative["converged"],
        "null_iterations": null["iterations"],
        "alternative_iterations": alternative["iterations"],
        "null_message": null["message"],
        "alternative_message": alternative["message"],
        "null_objective": null["objective"],
        "null_coefficients": null["parameters"][:-1].reshape(
            null_design.shape[1], n_paths - 1
        ),
        "alternative_coefficients": alternative["parameters"][
            : alternative_design.shape[1] * (n_paths - 1)
        ].reshape(alternative_design.shape[1], n_paths - 1),
    }
