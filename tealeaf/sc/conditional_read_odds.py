"""Experimental finite-count conditional local-read odds inference.

Within each subject/primer, condition a two-type, two-class table on margins.
Class capture ratios cancel only if shared by both cell types within a
subject/primer stratum. A Gaussian contrast shared by that
subject's primers is integrated by adaptive quadrature. The working random
slope distribution is conditional on margins, not an unconditional theorem
about RNA PSI. Validate with unconditional count nulls before use.
"""

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special, stats

MODEL_VERSION = "conditional_local_read_random_slope_v1"


@dataclass(frozen=True)
class ConditionalReadStratum:
    """Finite centered support and normalized central-hypergeometric mass."""

    centered_support: np.ndarray
    log_base: np.ndarray

    def evaluate(self, effects):
        """Log PMF of the observed table, score and information, Q-vectors."""
        effects = np.atleast_1d(np.asarray(effects, dtype=float))
        if effects.ndim != 1 or not np.isfinite(effects).all():
            raise ValueError("finite one-dimensional log-odds effects required")
        mass = self.log_base[None, :] + effects[:, None] * self.centered_support
        normalizer = special.logsumexp(mass, axis=1)
        probabilities = np.exp(mass - normalizer[:, None])
        average = probabilities @ self.centered_support
        information = np.sum(probabilities * (self.centered_support - average[:, None])**2, axis=1)
        observed_log_base = self.log_base[self.centered_support == 0]
        if len(observed_log_base) != 1:
            raise ValueError("observed table must be in its finite support")
        return observed_log_base[0] - normalizer, -average, information


class ConditionalReadOdds:
    """M by P by 2 by 2 integer counts, type then included/excluded class.

    Each informative subject has at least one primer with both types present
    and both classes represented across those types. Uninformative strata
    contribute exactly one to the conditional likelihood. Their subjects are
    recorded but cannot create an extra measured subject or test degree.
    """

    def __init__(self, counts, *, maximum_support=200_001):
        counts = np.asarray(counts)
        if counts.ndim != 4 or counts.shape[2:] != (2, 2) or not len(counts) or not counts.shape[1] or not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any() or (counts > 2**50).any() or maximum_support < 2:
            raise ValueError("nonnegative exact M by P by 2 by 2 integer counts required")
        self.counts = counts.astype(np.int64, copy=True)
        self.subject_strata, self.informative_subjects = [], []
        for index, subject in enumerate(self.counts):
            strata = []
            for table in subject:
                first, second = map(int, table.sum(axis=1))
                included = int(table[:, 0].sum())
                lower, upper = max(0, included - first), min(included, second)
                if upper <= lower:
                    continue
                if upper - lower + 1 > maximum_support:
                    raise ValueError("conditional support exceeds declared resource limit")
                support = np.arange(lower, upper + 1, dtype=float)
                # Fixed-margin constants do not depend on the effect. Explicit
                # normalization removes large log-binomial coordinate offsets.
                base = stats.hypergeom.logpmf(support, first + second, included, second)
                base -= special.logsumexp(base)
                strata.append(ConditionalReadStratum(support - int(table[1, 0]), base))
            if strata:
                self.informative_subjects.append(index)
                self.subject_strata.append(tuple(strata))
        self.informative_subjects = np.asarray(self.informative_subjects, dtype=int)

    @property
    def n_informative_subjects(self):
        return len(self.informative_subjects)

    def subject_evaluate(self, index, effects):
        result = [stratum.evaluate(effects) for stratum in self.subject_strata[index]]
        return tuple(np.sum([value[field] for value in result], axis=0) for field in range(3))

    def marginal_objective(self, mean, sd, *, nodes=21, gradient=True):
        """Integrate shared subject log odds with normalized adaptive GH.

        Posterior log densities are strictly concave for positive SD. Brackets
        for the unique mode follow from the finite sufficient-statistic range.
        The data-score identities give mean and log-SD derivatives without
        subtracting nearly equal prior scores at the zero-variance limit.
        Quadrature weights are not renormalized as an importance shortcut.
        """
        if not np.isfinite([mean, sd]).all() or sd < 0 or nodes < 3 or not self.n_informative_subjects:
            raise ValueError("finite mean, nonnegative SD, quadrature and measured subjects required")
        objective, derivative = 0., np.zeros(2)
        abscissae, weights = special.roots_hermitenorm(nodes)
        log_weights = np.log(weights / np.sqrt(2 * np.pi))
        for index, strata in enumerate(self.subject_strata):
            if sd == 0:
                value, score, _ = self.subject_evaluate(index, [mean])
                objective -= value[0]
                derivative[0] -= score[0]
                continue
            score_lower = -sum(stratum.centered_support.max() for stratum in strata)
            score_upper = -sum(stratum.centered_support.min() for stratum in strata)
            lower, upper = mean + sd**2 * score_lower, mean + sd**2 * score_upper
            def posterior_score(effect):
                return float(self.subject_evaluate(index, [effect])[1][0] - (effect - mean) / sd**2)
            mode = optimize.brentq(posterior_score, lower, upper, xtol=1e-11, rtol=1e-12)
            information = self.subject_evaluate(index, [mode])[2][0]
            scale = 1 / np.sqrt(information + 1 / sd**2)
            effects = mode + scale * abscissae
            values, scores, _ = self.subject_evaluate(index, effects)
            joint = log_weights + values + np.log(scale / sd) + .5 * abscissae**2 - .5 * ((effects - mean) / sd)**2
            integral = special.logsumexp(joint)
            posterior = np.exp(joint - integral)
            objective -= integral
            derivative[0] -= posterior @ scores
            derivative[1] -= posterior @ (scores * (effects - mean))
        if not np.isfinite(objective) or not np.isfinite(derivative).all():
            raise ValueError("conditional marginal likelihood exceeds finite range")
        return (float(objective), derivative) if gradient else float(objective)


def conditional_read_odds_test(likelihood, *, nodes=21, max_iter=100, quadrature_tolerance=1e-4):
    """Nested approximate LRT, refit heterogeneity under null and alternative.

    A native chi-square tail is diagnostic, not count-null or family-FDR
    certification. Compare positive SD and exact zero; failures stay p=1.
    Doubling quadrature checks both likelihoods and the likelihood ratio.
    Bounds and nonconvergence are exported, not silently called biological fit.
    """
    if likelihood.n_informative_subjects < 4:
        raise ValueError("four actually informative subject tables required")
    mean_bound, lower, upper = 20., -8., np.log(12.)
    options = dict(maxiter=max_iter, ftol=1e-10, gtol=1e-6, maxls=40)
    null_zero = likelihood.marginal_objective(0., 0., nodes=nodes, gradient=False)
    def null_objective(parameters):
        value, derivative = likelihood.marginal_objective(0., np.exp(parameters[0]), nodes=nodes)
        return value, derivative[1:]
    null_fits = [optimize.minimize(null_objective, [np.log(sd)], jac=True, method="L-BFGS-B", bounds=[(lower, upper)], options=options) for sd in (.2, 1., 3.)]
    finite = [fit for fit in null_fits if fit.success and np.isfinite(fit.fun)]
    null = min(finite, key=lambda fit: fit.fun) if finite else None
    null_sd = float(np.exp(null.x[0])) if null is not None and null.fun < null_zero else 0.
    null_value = min(null_zero, null.fun) if null is not None else null_zero
    def alternative_zero(parameters):
        value, derivative = likelihood.marginal_objective(parameters[0], 0., nodes=nodes)
        return value, derivative[:1]
    zero = optimize.minimize(alternative_zero, [0.], jac=True, method="L-BFGS-B", bounds=[(-mean_bound, mean_bound)], options=options)
    def alternative_objective(parameters):
        return likelihood.marginal_objective(parameters[0], np.exp(parameters[1]), nodes=nodes)
    initial_mean = float(zero.x[0]) if zero.success else 0.
    starts = [(0., max(null_sd, .2)), (initial_mean, .5), (initial_mean, 2.)]
    fits = [optimize.minimize(alternative_objective, [mean, np.log(sd)], jac=True, method="L-BFGS-B", bounds=[(-mean_bound, mean_bound), (lower, upper)], options=options) for mean, sd in starts]
    candidates = [(float(fit.fun), float(fit.x[0]), float(np.exp(fit.x[1])), fit) for fit in fits if fit.success and np.isfinite(fit.fun)]
    if zero.success and np.isfinite(zero.fun):
        candidates.append((float(zero.fun), float(zero.x[0]), 0., zero))
    if not candidates:
        raise ValueError("alternative conditional marginal optimization failed")
    alt_value, alt_mean, alt_sd, alt_fit = min(candidates, key=lambda fit: fit[0])
    null_fine = likelihood.marginal_objective(0., null_sd, nodes=2 * nodes, gradient=False)
    alternative_fine = likelihood.marginal_objective(alt_mean, alt_sd, nodes=2 * nodes, gradient=False)
    coarse, fine = 2 * (null_value - alt_value), 2 * (null_fine - alternative_fine)
    error = max(abs(null_value - null_fine), abs(alt_value - alternative_fine), abs(coarse - fine))
    boundary = abs(alt_mean) > mean_bound - 1e-3 or max(null_sd, alt_sd) > np.exp(upper) * (1 - 1e-3)
    converged = bool(null is not None and not boundary and min(coarse, fine) >= -1e-6 and error <= quadrature_tolerance)
    statistic = max(float(fine), 0.) if converged else 0.
    return dict(model_version=MODEL_VERSION, p_value=float(stats.chi2.sf(statistic, 1)) if converged else 1., statistic=statistic, converged=converged, n_subjects=likelihood.n_informative_subjects, n_requested_subjects=len(likelihood.counts), log_odds_effect=alt_mean if converged else np.nan, null_subject_sd=null_sd, alternative_subject_sd=alt_sd, null_objective=float(null_fine), alternative_objective=float(alternative_fine), quadrature_error=float(error), parameter_boundary=bool(boundary), null_optimizer_success=null is not None, alternative_optimizer_success=bool(alt_fit.success), scope="native-tail conditional local-read odds diagnostic, not RNA PSI or certified count-null inference")


def simulate_local_read_counts(totals, baseline_logits, effects, rng):
    """Unconditional binomial counts, M by P by two types by two classes.

    M-vector log-odds effects are shared across primers. Baselines M by P
    absorb path usage and primer capture. Totals M by P by 2 remain fixed.
    This is a local-marker working law, not an EC/read alignment simulator.
    """
    totals, baseline, effects = np.asarray(totals), np.asarray(baseline_logits, dtype=float), np.asarray(effects, dtype=float)
    if totals.ndim != 3 or totals.shape[2] != 2 or not len(totals) or baseline.shape != totals.shape[:2] or effects.shape != (len(totals),) or not np.isfinite(totals).all() or (totals < 0).any() or (totals != np.floor(totals)).any() or (totals > 2**50).any() or not np.isfinite(baseline).all() or not np.isfinite(effects).all():
        raise ValueError("integer type totals and aligned finite baseline/effects required")
    totals = totals.astype(np.int64)
    probabilities = special.expit(baseline[:, :, None] + effects[:, None, None] * np.array([-.5, .5]))
    included = rng.binomial(totals, probabilities)
    return np.stack([included, totals - included], axis=-1)


def simulate_conditional_read_counts(counts, effects, rng):
    """Draw noncentral hypergeometric tables, preserving every observed margin.

    This is the conditional working law, distinct from generating independent
    binomials with random effects and then conditioning their observations.
    Effects must be assigned independently of these fixed margins.
    """
    original = ConditionalReadOdds(counts).counts
    effects = np.asarray(effects, dtype=float)
    if effects.shape != (len(original),) or not np.isfinite(effects).all() or (np.abs(effects) > 700).any():
        raise ValueError("one finite supported log-odds effect per requested subject required")
    generated = original.copy()
    for subject, local in enumerate(original):
        for primer, table in enumerate(local):
            first, second = map(int, table.sum(axis=1))
            included = int(table[:, 0].sum())
            lower, upper = max(0, included - first), min(included, second)
            if upper <= lower:
                selected = lower
            elif effects[subject] == 0:
                selected = int(rng.hypergeometric(included, first + second - included, second))
            else:
                selected = int(stats.nchypergeom_fisher.rvs(first + second, included, second, np.exp(effects[subject]), random_state=rng))
            generated[subject, primer] = [[included - selected, first - included + selected], [selected, second - selected]]
    return generated
