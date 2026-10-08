"""Experimental unconditional binomial local-read random-intercept/slope model.

Primer-specific fixed baselines and a shared subject baseline account for
capture differences that are shared by cell types. Independent Gaussian
subject slopes model heterogeneous local read odds. This is not absolute RNA
PSI and does not cancel type-specific capture bias or unmodeled marker paths.
"""

import numpy as np
from scipy import special, optimize, stats

MODEL_VERSION = "local_read_binomial_random_intercept_slope_v1"


class LocalReadMixed:
    """M by P by two types by included/excluded integer counts, no pseudocount.

    Latent standard normals z_u have dimension two. Predictors are
    gamma_p + sd_baseline*z_u0 + s_c*(effect + sd_slope*z_u1)/2,
    with s_c=-1,+1. Means have no penalty. Subject baseline and slope are
    independent in this experimental model; test correlated-null sensitivity
    separately before real adoption. Zero-observation strata are not evidence.
    """

    def __init__(self, counts):
        counts = np.asarray(counts)
        if counts.ndim != 4 or counts.shape[2:] != (2, 2) or not len(counts) or not counts.shape[1] or not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any() or (counts > 2**50).any():
            raise ValueError("nonnegative exact M by P by two by two counts required")
        self.counts = counts.astype(np.int64, copy=True)
        self.active_primers = np.flatnonzero(self.counts.sum(axis=(0, 2, 3)) > 0)
        if not len(self.active_primers):
            raise ValueError("local read evidence required")
        kept = self.counts[:, self.active_primers]
        self.included = kept[..., 0].reshape(len(kept), -1)
        self.excluded = kept[..., 1].reshape(len(kept), -1)
        self.totals = self.included + self.excluded
        self.primers = len(self.active_primers)
        self.signs = np.tile(np.array([-.5, .5]), self.primers)
        self.primer_design = np.repeat(np.eye(self.primers), 2, axis=0)
        design = np.column_stack([self.primer_design, self.signs])
        if np.linalg.matrix_rank(design[self.totals.sum(axis=0) > 0]) < self.primers + 1:
            raise ValueError("primer and cell-type means are confounded")
        self.n_paired_subjects = int((kept.sum(axis=(1, 3)) > 0).all(axis=1).sum())

    def _data(self, predictors):
        positive, negative = special.expit(predictors), special.expit(-predictors)
        included, excluded = self.included[:, None, :], self.excluded[:, None, :]
        nll = np.sum(included * np.logaddexp(0., -predictors) + excluded * np.logaddexp(0., predictors), axis=-1)
        residual = included * negative - excluded * positive
        information = (included + excluded) * positive * negative
        return nll, residual, information

    def marginal_objective(self, primer_means, effect, baseline_sd, slope_sd, *, nodes=9, gradient=True):
        """Batched two-dimensional adaptive Gaussian quadrature over subjects.

        Posterior modes use standard-normal coordinates, so the zero-variance
        cases remain well conditioned and are evaluated exactly in that latent
        direction. Returns objective up to fixed binomial constants and a
        P+3 gradient, primer means, effect, log baseline SD, log slope SD.
        """
        means = np.asarray(primer_means, dtype=float)
        if means.shape != (self.primers,) or not np.isfinite(means).all() or not np.isfinite([effect, baseline_sd, slope_sd]).all() or min(baseline_sd, slope_sd) < 0 or nodes < 3:
            raise ValueError("aligned finite means and nonnegative latent SDs required")
        fixed = self.primer_design @ means + self.signs * effect
        latent = np.column_stack([np.full(2 * self.primers, baseline_sd), self.signs * slope_sd])
        modes = np.zeros((len(self.counts), 2))
        identity = np.eye(2)
        def mode_state(values):
            nll, residual, information = self._data((fixed[None, :] + values @ latent.T)[:, None, :])
            objective = nll[:, 0] + .5 * np.sum(values**2, axis=1)
            derivative = values - residual[:, 0] @ latent
            curvature = identity[None, :, :] + np.einsum("mr,ri,rj->mij", information[:, 0], latent, latent)
            return objective, derivative, curvature
        for _ in range(80):
            objective, derivative, curvature = mode_state(modes)
            step = np.linalg.solve(curvature, derivative[..., None])[..., 0]
            if np.max(np.abs(step)) < 1e-10:
                break
            scales = np.ones(len(modes))
            decrement = np.sum(derivative * step, axis=1)
            for _ in range(40):
                candidates = modes - scales[:, None] * step
                proposal = mode_state(candidates)[0]
                failed = proposal > objective - 1e-4 * scales * decrement + 1e-12
                if not failed.any():
                    break
                scales[failed] *= .5
            else:
                raise ValueError("subject posterior Newton line search failed")
            modes = candidates
        else:
            raise ValueError("subject posterior modes did not converge")
        _, _, curvature = mode_state(modes)
        cholesky = np.linalg.cholesky(curvature)
        proposal_factor = np.linalg.inv(cholesky).transpose(0, 2, 1)
        abscissae, weights = special.roots_hermitenorm(nodes)
        grid = np.stack(np.meshgrid(abscissae, abscissae, indexing="ij"), axis=-1).reshape(-1, 2)
        log_weights = (np.log(weights)[:, None] + np.log(weights)[None, :] - np.log(2 * np.pi)).reshape(-1)
        positions = modes[:, None, :] + np.einsum("qj,mkj->mqk", grid, proposal_factor)
        predictors = fixed[None, None, :] + positions @ latent.T
        nll, residual, _ = self._data(predictors)
        log_determinant = -np.log(np.diagonal(cholesky, axis1=1, axis2=2)).sum(axis=1)
        joint = log_weights[None, :] + log_determinant[:, None] + .5 * np.sum(grid**2, axis=1)[None, :] - .5 * np.sum(positions**2, axis=2) - nll
        evidence = special.logsumexp(joint, axis=1)
        objective = -evidence.sum()
        if not gradient:
            return float(objective)
        posterior = np.exp(joint - evidence[:, None])
        weighted = posterior[:, :, None] * residual
        derivative = np.r_[np.sum(weighted, axis=(0, 1)) @ self.primer_design, np.sum(weighted * self.signs[None, None, :]), baseline_sd * np.sum(weighted * positions[:, :, 0, None]), slope_sd * np.sum(weighted * positions[:, :, 1, None] * self.signs[None, None, :])]
        if not np.isfinite(objective) or not np.isfinite(derivative).all():
            raise ValueError("binomial marginal objective exceeds finite range")
        return float(objective), -derivative


def local_read_mixed_test(likelihood, *, nodes=9, max_iter=150, quadrature_tolerance=1e-3):
    """Native chi-square LRT with both latent variances refitted and checked.

    Compare all four exact zero/positive variance regimes. Doubled-order
    likelihood validation is required. Native tails still require actual
    count-null assessment and are not certified biological p-values.
    """
    if likelihood.n_paired_subjects < 4:
        raise ValueError("four subjects with local reads in both types required")
    counts = likelihood.counts[:, likelihood.active_primers]
    pooled = counts.sum(axis=(0, 2))
    initial = special.logit((pooled[:, 0] + .5) / (pooled.sum(axis=1) + 1))
    options = dict(maxiter=max_iter, maxls=40, ftol=1e-10, gtol=1e-6)
    def fit_family(alternative, starts=None):
        results = []
        for flags in ((False, False), (True, False), (False, True), (True, True)):
            free = likelihood.primers + int(alternative)
            means = np.r_[initial, 0.] if alternative else initial.copy()
            scales = [.7, .5]
            if starts is not None:
                means = np.r_[starts[1], 0.] if alternative else starts[1].copy()
                scales = [max(starts[3], .2), max(starts[4], .2)]
            parameters = np.r_[means, [np.log(scales[index]) for index, active in enumerate(flags) if active]]
            def unpack(values):
                primer = values[:likelihood.primers]
                effect = values[likelihood.primers] if alternative else 0.
                sd, position = [0., 0.], free
                for index, active in enumerate(flags):
                    if active:
                        sd[index] = np.exp(values[position])
                        position += 1
                return primer, effect, *sd
            def objective(values):
                primer, effect, first, second = unpack(values)
                value, derivative = likelihood.marginal_objective(primer, effect, first, second, nodes=nodes)
                selected = list(range(likelihood.primers)) + ([likelihood.primers] if alternative else []) + [likelihood.primers + 1 + index for index, active in enumerate(flags) if active]
                return value / len(likelihood.counts), derivative[selected] / len(likelihood.counts)
            try:
                fit = optimize.minimize(objective, parameters, jac=True, method="L-BFGS-B", bounds=[(-20., 20.)] * free + [(-8., np.log(12.))] * sum(flags), options=options)
            except (ValueError, np.linalg.LinAlgError):
                continue
            if fit.success and np.isfinite(fit.fun):
                primer, effect, first, second = unpack(fit.x)
                results.append((float(fit.fun * len(likelihood.counts)), primer, effect, first, second, fit))
        if not results:
            raise ValueError("binomial random-intercept/slope fitting failed")
        return min(results, key=lambda result: result[0])
    null = fit_family(False)
    alternative = fit_family(True, starts=null)
    null_fine = likelihood.marginal_objective(*null[1:5], nodes=2 * nodes, gradient=False)
    alternative_fine = likelihood.marginal_objective(*alternative[1:5], nodes=2 * nodes, gradient=False)
    coarse, fine = 2 * (null[0] - alternative[0]), 2 * (null_fine - alternative_fine)
    error = max(abs(null[0] - null_fine), abs(alternative[0] - alternative_fine), abs(coarse - fine))
    boundary = max(np.max(np.abs(null[1])), np.max(np.abs(alternative[1])), abs(alternative[2])) > 20 - 1e-3 or max(null[3], null[4], alternative[3], alternative[4]) > 12 * (1 - 1e-3)
    converged = bool(not boundary and error <= quadrature_tolerance and min(coarse, fine) >= -1e-6)
    statistic = max(float(fine), 0.) if converged else 0.
    return dict(model_version=MODEL_VERSION, p_value=float(stats.chi2.sf(statistic, 1)) if converged else 1., statistic=statistic, converged=converged, n_subjects=likelihood.n_paired_subjects, n_requested_subjects=len(likelihood.counts), log_odds_effect=float(alternative[2]) if converged else np.nan, null_baseline_sd=float(null[3]), null_subject_sd=float(null[4]), alternative_baseline_sd=float(alternative[3]), alternative_subject_sd=float(alternative[4]), null_objective=float(null_fine), alternative_objective=float(alternative_fine), quadrature_error=float(error), parameter_boundary=bool(boundary), scope="native-tail unconditional binomial local-read odds model with independent normal baseline/slope, not RNA PSI or certified biological inference")
