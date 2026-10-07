"""Experimental paired EC scores allowing unknown shared read opportunities.

Condition on each primer's EC total across the two types. Multiplicative EC
opportunity factors shared by those types cancel. Type-specific opportunity
bias does not cancel; this is not an absolute-usage quantification model.
"""

import numpy as np
from scipy import linalg, optimize, sparse, special

from . import differential
from .path_bias import SharedPathNullProblem, SharedPathNullFit

MODEL_VERSION = "conditional_ec_opportunity_v1"


def binary_fragment_opportunity_kernels(designs):
    """Undo EC-category degree normalization under the paired binary model.

    Requires global (not gene-sliced) oligoDT-TPM/random-hexamer-length maps,
    with identical binary support and column-constant positive entries. Let
    d_t be global EC degree and L_t the hexamer column sum, its relative
    effective length. Multiplying both maps by d_t/L_t gives membership/L_t
    for oligoDT and membership for hexamers. Shared EC row opportunities are
    then left unknown. This assumes uniform read starts conditional on primer
    sampling, not arbitrary transcript-specific positional bias.
    """
    if len(designs) != 2:
        raise ValueError("two global primer maps required")
    poly, hexamer = [sparse.csr_matrix(value, dtype=float, copy=True) for value in designs]
    for mapping in (poly, hexamer):
        mapping.eliminate_zeros()
        mapping.sort_indices()
        if not np.isfinite(mapping.data).all() or (mapping.data < 0).any():
            raise ValueError("finite nonnegative binary maps required")
    if poly.shape != hexamer.shape or not np.array_equal(poly.indptr, hexamer.indptr) or not np.array_equal(poly.indices, hexamer.indices):
        raise ValueError("binary primer maps must have identical transcript support")
    degree = np.asarray(poly.getnnz(axis=0)).ravel()
    poly_sums = np.asarray(poly.sum(axis=0)).ravel()
    lengths = np.asarray(hexamer.sum(axis=0)).ravel()
    supported = degree > 0
    if not np.allclose(poly_sums[supported], 1., atol=1e-8) or (lengths[supported] <= 0).any():
        raise ValueError("global oligoDT-TPM column normalization required")
    for mapping, sums in ((poly, poly_sums), (hexamer, lengths)):
        expected = sums[mapping.indices] / degree[mapping.indices]
        if not np.allclose(mapping.data, expected, rtol=1e-8, atol=1e-12):
            raise ValueError("column-constant binary likelihoods required, weighted designs are not supported")
    factors = np.ones(poly.shape[1])
    factors[supported] = degree[supported] / lengths[supported]
    return tuple((mapping @ sparse.diags(factors)).tocsr() for mapping in (poly, hexamer))


class ConditionalPathNullProblem:
    """Shared S-path null, two free transcript mixtures and primer intercepts.

    Counts are 2 by K_p and maps K_p by T. Parameters contain the original
    (S-1)+2*(T-S) composition coordinates and P active-primer log-depth ratios.
    Null path proportions are not necessarily identifiable from conditional
    counts. Count-null and observed-power checks are required before use.
    """

    def __init__(self, counts, designs, baseline, path_index, nuisance_pseudocount=1e-4):
        self.structure = SharedPathNullProblem(counts, designs, baseline, path_index, nuisance_pseudocount)
        if self.structure.types != 2:
            raise ValueError("conditional event prototype requires exactly two types")
        self.active_primers = [index for index, totals in enumerate(self.structure.totals) if (totals > 0).all()]
        if not self.active_primers:
            raise ValueError("conditional contrast needs a primer observed in both types")
        offsets = [np.log(self.structure.totals[index][1] / self.structure.totals[index][0]) for index in self.active_primers]
        self.initial = np.r_[self.structure.initial, offsets]
        self.dimension = len(self.initial)

    def objective(self, parameters):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (self.dimension,):
            raise ValueError("conditional null parameter dimension differs")
        theta, jacobian, _, value, prior_gradient = self.structure.composition(parameters[:self.structure.dimension])
        gradient = np.r_[prior_gradient, np.zeros(len(self.active_primers))]
        for offset, primer in enumerate(self.active_primers):
            counts, mapping = self.structure.counts[primer], self.structure.maps[primer]
            trials = counts.sum(axis=0)
            active = trials > 0
            mass = mapping[active] @ theta.T
            if (mass <= 0).any():
                return np.inf, np.zeros(self.dimension)
            predictor = parameters[self.structure.dimension + offset] + np.log(mass[:, 1]) - np.log(mass[:, 0])
            probability = special.expit(predictor)
            observed = counts[1, active]
            value += float(trials[active] @ np.logaddexp(0., predictor) - observed @ predictor)
            residual = trials[active] * probability - observed
            derivative = ((mapping[active] * theta[1]) @ jacobian[1]) / mass[:, 1, None] - ((mapping[active] * theta[0]) @ jacobian[0]) / mass[:, 0, None]
            gradient[:self.structure.dimension] += derivative.T @ residual
            gradient[self.structure.dimension + offset] += residual.sum()
        return value, gradient

    def fit(self, *, max_iter=300, tolerance=1e-12, multistart=False):
        scale = float(sum(self.structure.counts[index].sum() for index in self.active_primers))

        def objective(parameters):
            value, gradient = self.objective(parameters)
            return value / scale, gradient / scale

        starts = [self.initial]
        if multistart:
            starts.append(np.r_[np.zeros(self.structure.dimension), self.initial[self.structure.dimension:]])
        results = [optimize.minimize(objective, initial, jac=True, method="L-BFGS-B", bounds=[(-30., 30.)] * self.dimension, options=dict(maxiter=int(max_iter), maxls=50, ftol=float(tolerance), gtol=1e-8)) for initial in starts]
        finite = [index for index, result in enumerate(results) if np.isfinite(result.fun)]
        selected = min(finite, key=lambda index: results[index].fun) if finite else 0
        if finite:
            best = results[selected].fun
            tied = [index for index in finite if results[index].fun <= best + tolerance * max(1., abs(best))]
            selected = min(tied, key=lambda index: (not results[index].success, results[index].fun))
        result = results[selected]
        theta, _, psi, _, _ = self.structure.composition(result.x[:self.structure.dimension])
        fit = SharedPathNullFit(psi, theta, bool(result.success), float(result.fun * scale), int(result.nit), float(linalg.norm(result.jac, ord=np.inf) * scale), str(result.message), len(starts), selected)
        return fit, result.x[self.structure.dimension:]


def efficient_conditional_path_score(problem, theta, offsets):
    """Efficient S-1 ILR score, information, Dirichlet shape and reference.

    Nuisance projection uses linear transcript-simplex directions, a common
    path coordinate and a separate intercept per active primer. Unknown
    shared row opportunity factors cancel from every predictor derivative.
    """
    structure = problem.structure
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (2, len(structure.baseline)) or np.asarray(offsets).shape != (len(problem.active_primers),):
        raise ValueError("aligned conditional composition and primer offsets required")
    paths, size, basis = structure.path_index, structure.size, structure.basis
    psi = differential.path_proportions(theta[0], paths)
    if (psi <= 0).any() or not np.allclose(differential.path_proportions(theta[1], paths), psi, atol=1e-8):
        raise ValueError("positive shared path proportions required")
    columns = []
    groups = [np.flatnonzero(paths == path) for path in range(size)]
    outside = np.flatnonzero(paths < 0)
    for group in groups + ([outside] if len(outside) else []):
        if len(group) > 1:
            embedded = np.zeros((len(paths), len(group) - 1))
            embedded[group] = differential.helmert_basis(len(group))
            columns.extend(embedded.T)
    if len(outside):
        direction = np.zeros(len(paths))
        for path, group in enumerate(groups):
            direction[group] = psi[path] / len(group)
        direction[outside] = -1 / len(outside)
        columns.append(direction / linalg.norm(direction))
    nuisance = np.asarray(columns).T if columns else np.zeros((len(paths), 0))
    target_jacobian = np.zeros((len(paths), size - 1))
    target_jacobian[paths >= 0] = basis[paths[paths >= 0]] - psi @ basis
    blocks, nulls, residuals = [], [], []
    for offset, primer in enumerate(problem.active_primers):
        counts, mapping = structure.counts[primer], structure.maps[primer]
        trials = counts.sum(axis=0)
        active = trials > 0
        mapping, trials = mapping[active], trials[active]
        mass = mapping @ theta.T
        if (mass <= 0).any():
            raise ValueError("positive conditional counts require positive mapped mass")
        predictor = offsets[offset] + np.log(mass[:, 1]) - np.log(mass[:, 0])
        probability = special.expit(predictor)
        variance = trials * probability * special.expit(-predictor)
        if (variance <= 0).any():
            raise ValueError("conditional binomial information is zero")
        target_a = (mapping @ (theta[0, :, None] * target_jacobian)) / mass[:, 0, None]
        target_b = (mapping @ (theta[1, :, None] * target_jacobian)) / mass[:, 1, None]
        intercepts = np.zeros((len(trials), len(problem.active_primers)))
        intercepts[:, offset] = 1.
        null = np.column_stack((target_b - target_a, -(mapping @ nuisance) / mass[:, 0, None], (mapping @ nuisance) / mass[:, 1, None], intercepts))
        blocks.append(np.sqrt(variance)[:, None] * target_b)
        nulls.append(np.sqrt(variance)[:, None] * null)
        residuals.append((counts[1, active] - trials * probability) / np.sqrt(variance))
    target, null, residual = np.vstack(blocks), np.vstack(nulls), np.concatenate(residuals)
    norms = linalg.norm(null, axis=0)
    nonzero = norms > np.finfo(float).tiny
    orthogonal = linalg.orth(null[:, nonzero] / norms[nonzero], rcond=1e-10)
    efficient = target - orthogonal @ (orthogonal.T @ target)
    score, information = efficient.T @ residual, efficient.T @ efficient
    shape = 2 * (basis.T / psi) @ basis
    return score, (information + information.T) / 2, shape, target.T @ target
