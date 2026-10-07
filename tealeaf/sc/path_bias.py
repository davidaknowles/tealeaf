"""Experimental removal of coverage-dependent path-smoothing contrasts.

The correction fits a common target-path null but permits independent
within-path transcript mixtures and outside-block mass at each type. It then
subtracts the SAME smoothed estimator evaluated on null expected EC counts.
This deterministic correction is not an exact finite-count bias estimate;
count-null calibration is required before interpreting any downstream test.
Production inference and reporting defaults are unchanged.
"""

from dataclasses import dataclass

import numpy as np
from scipy import optimize, sparse, special

from . import differential
from .ec_block_glmm import pooled_isoform_weights


@dataclass
class SharedPathNullFit:
    """Common S-path proportions and C by T nuisance-specific mixtures."""

    path_proportions: np.ndarray
    theta: np.ndarray
    converged: bool
    objective: float
    iterations: int
    gradient_norm: float


class SharedPathNullProblem:
    """Primer-conditioned null with (S-1)+C*(T-S) free coordinates.

    Primer count matrices are C by K_p; maps are K_p by T. Path indices
    have length T, consecutive 0..S-1 inside the block and -1 outside.
    The shared path proportions have no prior. A tiny total pseudocount
    stabilizes each within-path/outside mixture, never the target contrast.
    """

    def __init__(self, counts, designs, baseline, path_index, nuisance_pseudocount=1e-4):
        self.baseline = np.asarray(baseline, dtype=float)
        self.path_index = np.asarray(path_index, dtype=int)
        if self.baseline.ndim != 1 or self.path_index.shape != self.baseline.shape or not np.isfinite(self.baseline).all() or (self.baseline < 0).any() or self.baseline.sum() <= 0 or (self.path_index < -1).any():
            raise ValueError("finite nonnegative baseline and aligned path indices required")
        self.baseline = np.maximum(self.baseline, 1e-12)
        self.baseline /= self.baseline.sum()
        paths = np.unique(self.path_index[self.path_index >= 0])
        if len(paths) < 2 or not np.array_equal(paths, np.arange(len(paths))):
            raise ValueError("at least two consecutive local paths required")
        self.size = len(paths)
        self.basis = differential.helmert_basis(self.size)
        self.counts = tuple(np.asarray(value, dtype=float) for value in counts)
        self.maps = tuple(np.asarray(sparse.csr_matrix(value).todense()) for value in designs)
        if not self.counts or len(self.counts) != len(self.maps) or self.counts[0].ndim != 2:
            raise ValueError("aligned primer count matrices and mappings required")
        self.types = self.counts[0].shape[0]
        if self.types < 2 or any(value.shape != (self.types, mapping.shape[0]) or mapping.shape[1] != len(self.baseline) or not np.isfinite(value).all() or not np.isfinite(mapping).all() or (value < 0).any() or (mapping < 0).any() for value, mapping in zip(self.counts, self.maps)):
            raise ValueError("finite nonnegative C by K counts and K by T maps required")
        self.totals = tuple(value.sum(axis=1) for value in self.counts)
        if (sum(self.totals) <= 0).any() or not np.isfinite(nuisance_pseudocount) or nuisance_pseudocount < 0:
            raise ValueError("positive total per type and nonnegative nuisance prior required")
        self.nuisance_pseudocount = float(nuisance_pseudocount)
        self.groups = [np.flatnonzero(self.path_index == path) for path in paths]
        self.outside = np.flatnonzero(self.path_index < 0)
        if len(self.outside):
            self.groups.append(self.outside)
        self.group_bases = [differential.helmert_basis(len(group)) if len(group) > 1 else np.zeros((1, 0)) for group in self.groups]
        self.nuisance_dimension = len(self.baseline) - self.size
        self.dimension = self.size - 1 + self.types * self.nuisance_dimension
        psi = differential.path_proportions(self.baseline, self.path_index)
        nuisance = [basis.T @ np.log(self.baseline[group] / self.baseline[group].sum()) for group, basis in zip(self.groups, self.group_bases)]
        if len(self.outside):
            h = self.baseline[self.path_index >= 0].sum()
            nuisance.append(np.array([special.logit(h)]))
        nuisance = np.concatenate(nuisance)
        self.initial = np.r_[self.basis.T @ np.log(psi), np.tile(nuisance, self.types)]

    def composition(self, parameters):
        """Return theta C by T, log-theta derivative C by T by P, psi S."""
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (self.dimension,):
            raise ValueError("parameter dimension differs from shared-path null")
        psi = special.softmax(self.basis @ parameters[:self.size - 1])
        theta = np.zeros((self.types, len(self.baseline)))
        jacobian = np.zeros((*theta.shape, self.dimension))
        penalty, penalty_gradient = 0., np.zeros(self.dimension)
        for cell in range(self.types):
            start = self.size - 1 + cell * self.nuisance_dimension
            offset = start
            h = special.expit(parameters[start + self.nuisance_dimension - 1]) if len(self.outside) else 1.
            for path, (group, basis) in enumerate(zip(self.groups, self.group_bases)):
                dimension = basis.shape[1]
                shares = special.softmax(basis @ parameters[offset:offset + dimension])
                theta[cell, group] = shares * (h * psi[path] if path < self.size else 1 - h)
                jacobian[cell, group, offset:offset + dimension] = basis - shares @ basis
                if path < self.size:
                    jacobian[cell, group, :self.size - 1] = self.basis[path] - psi @ self.basis
                if len(self.outside):
                    jacobian[cell, group, start + self.nuisance_dimension - 1] = 1 - h if path < self.size else -h
                if dimension:
                    penalty -= self.nuisance_pseudocount * float(np.log(shares).mean())
                    penalty_gradient[offset:offset + dimension] += self.nuisance_pseudocount * (shares @ basis)
                offset += dimension
        return theta, jacobian, psi, penalty, penalty_gradient

    def objective(self, parameters):
        theta, jacobian, _, value, gradient = self.composition(parameters)
        for counts, mapping, totals in zip(self.counts, self.maps, self.totals):
            for cell, total in enumerate(totals):
                if total <= 0:
                    continue
                mass = mapping @ theta[cell]
                normalizer = mass.sum()
                if normalizer <= 0 or ((counts[cell] > 0) & (mass <= 0)).any():
                    return np.inf, np.zeros(self.dimension)
                safe = np.maximum(mass, 1e-300)
                value -= float(counts[cell] @ np.log(safe)) - total * np.log(normalizer)
                score = theta[cell] * (-mapping.T @ (counts[cell] / safe) + total * mapping.sum(axis=0) / normalizer)
                gradient += jacobian[cell].T @ score
        return value, gradient

    def fit(self, *, max_iter=300, tolerance=1e-12):
        result = optimize.minimize(self.objective, self.initial, jac=True, method="L-BFGS-B", bounds=[(-30., 30.)] * self.dimension, options={"maxiter": int(max_iter), "ftol": float(tolerance), "gtol": 1e-8})
        theta, _, psi, _, _ = self.composition(result.x)
        return SharedPathNullFit(psi, theta, bool(result.success), float(result.fun), int(result.nit), float(np.linalg.norm(result.jac, ord=np.inf)))


def expected_ec_counts(counts, designs, theta):
    """Null expectations C by K_p, preserving every type/primer total."""
    result = []
    for observed, mapping in zip(counts, designs):
        observed = np.asarray(observed, dtype=float)
        mass = np.asarray(sparse.csr_matrix(mapping) @ np.asarray(theta).T).T
        totals, normalizers = observed.sum(axis=1), mass.sum(axis=1)
        if ((totals > 0) & (normalizers <= 0)).any():
            raise ValueError("observed primer has no expected model mass")
        result.append(mass * np.divide(totals, normalizers, out=np.zeros_like(totals), where=normalizers > 0)[:, None])
    return tuple(result)


def null_corrected_path_responses(data, path_index, labels, subjects, *, baseline=None, concentration=32., max_iter=300, reporting_concentration=None):
    """Return smoothed-minus-null responses for general C-type designs.

    Values are N_aggregate by (S-1) linear-proportion contrasts, NOT ILRs
    or simplex-constrained reported usage estimates. A fit failure invalidates
    the entire hypothesis rather than silently selecting successful subjects.
    Missing types and zero-total aggregates are structural exclusions only.
    """
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    if labels.shape != subjects.shape or labels.shape != (len(data.counts[0]),) or not np.isfinite(concentration) or concentration < 0:
        raise ValueError("aligned labels/subjects and finite nonnegative concentration required")
    if reporting_concentration is not None and (not np.isfinite(reporting_concentration) or reporting_concentration < 0):
        raise ValueError("finite nonnegative reporting concentration required")
    levels, encoded = np.unique(labels, return_inverse=True)
    if len(levels) < 2:
        raise ValueError("at least two types required")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    size = len(np.unique(np.asarray(path_index)[np.asarray(path_index) >= 0]))
    basis = differential.helmert_basis(size)
    values, retained_subjects, retained_labels, observed_psi, null_psi, null_fits, reporting_psi = [], [], [], [], [], [], []
    for subject in np.unique(subjects):
        local_levels = np.unique(encoded[subjects == subject])
        counts = [np.asarray([np.asarray(matrix[(subjects == subject) & (encoded == level)], dtype=float).sum(axis=0) for level in local_levels]) for matrix in data.counts]
        positive = sum(value.sum(axis=1) for value in counts) > 0
        local_levels, counts = local_levels[positive], tuple(value[positive] for value in counts)
        if len(local_levels) < 2:
            continue
        problem = SharedPathNullProblem(counts, data.compatibility, baseline, path_index)
        shared = problem.fit(max_iter=max_iter)
        if not shared.converged:
            raise ValueError(f"shared-path null optimization failed for subject {subject}")
        expected = expected_ec_counts(counts, data.compatibility, shared.theta)
        for local, level in enumerate(local_levels):
            fits = [differential.fit_free_isoform_paths(tuple(value[local] for value in source), data.compatibility, baseline, path_index, path_pseudocount=concentration, path_pseudocount_scaling="total", max_iter=max_iter, tolerance=1e-12) for source in (counts, expected)]
            if not all(fit.converged for fit in fits):
                raise ValueError(f"observed/null-expectation optimization failed for subject {subject}, type {levels[level]}")
            observed, reference = (fit.path_proportions for fit in fits)
            residual = observed - reference
            # Deterministically exact expectations should not produce a
            # biological contrast from sub-nanounit optimizer roundoff.
            residual[np.abs(residual) < 1e-9] = 0.
            values.append(basis.T @ residual)
            observed_psi.append(observed)
            null_psi.append(reference)
            if reporting_concentration is None or reporting_concentration == concentration:
                reporting_psi.append(observed)
            else:
                report = differential.fit_free_isoform_paths(tuple(value[local] for value in counts), data.compatibility, baseline, path_index, path_pseudocount=reporting_concentration, path_pseudocount_scaling="total", max_iter=max_iter, tolerance=1e-12)
                # Reporting failures do not alter successful inference.
                reporting_psi.append(report.path_proportions if report.converged else np.full(size, np.nan))
            retained_subjects.append(subject)
            retained_labels.append(level)
        null_fits.append(shared)
    return {"values": np.asarray(values).reshape(-1, size - 1), "subjects": np.asarray(retained_subjects), "encoded_labels": np.asarray(retained_labels), "levels": tuple(levels), "observed_proportions": np.asarray(observed_psi).reshape(-1, size), "null_smoothed_proportions": np.asarray(null_psi).reshape(-1, size), "reporting_proportions": np.asarray(reporting_psi).reshape(-1, size), "null_fits": null_fits, "concentration": concentration, "reporting_concentration": reporting_concentration}


def paired_null_corrected_path_test(data, path_index, labels, subjects, **kwargs):
    """Experimental paired t/Hotelling test on null-corrected proportions."""
    if len(np.unique(labels)) != 2:
        raise ValueError("paired null correction requires exactly two types")
    responses = null_corrected_path_responses(data, path_index, labels, subjects, **kwargs)
    subject_ids = np.unique(responses["subjects"])
    differences = np.asarray([responses["values"][(responses["subjects"] == subject) & (responses["encoded_labels"] == 1)][0] - responses["values"][(responses["subjects"] == subject) & (responses["encoded_labels"] == 0)][0] for subject in subject_ids]).reshape(-1, responses["values"].shape[1])
    tested = differential.paired_mean_test(differences)
    if len(subject_ids) < 4:
        tested.update(p_value=1., statistic=0., converged=False)
    return {**tested, **responses, "differences": differences, "subject_ids": subject_ids}
