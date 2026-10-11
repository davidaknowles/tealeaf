"""Compiled kernels for the grid random-subject binary path model.

Numba versions of ``path_marginal_grid.row_integral``, the per-subject
posterior-mode search and the adaptive Gauss-Hermite objective, plus the
profiled-precursor row likelihood of ``local_path_test``. They implement the
same formulas as the numpy references, which the tests compare against; the
gain is removing per-row and per-subject Python overhead.

Grid rows are packed into flat arrays: for row n, points start[n]:stop[n] of
log_p, log_q (log psi, log(1 - psi) at each grid t) and base (shifted row
log-likelihood plus trapezoid log weight), with first[n] and last[n] the
shifted log-likelihood at the grid ends for the analytic tails. Rows with
kind[n] = 1 use the exact Beta-binomial terms[n] = (path-0 count, path-1
count, log constant).
"""

import math

import numpy as np
from numba import njit

LIMIT = 18.


@njit(cache=True, error_model="numpy")
def digamma(x):
    result = 0.
    while x < 10.:
        result -= 1. / x
        x += 1.
    f = 1. / (x * x)
    return result + math.log(x) - .5 / x - f * (1. / 12 - f * (1. / 120 - f * (1. / 252 - f * (1. / 240 - f / 132))))


@njit(cache=True, error_model="numpy")
def trigamma(x):
    result = 0.
    while x < 10.:
        result += 1. / (x * x)
        x += 1.
    f = 1. / (x * x)
    return result + 1. / x + .5 * f + f / x * (1. / 6 - f * (1. / 30 - f * (1. / 42 - f * (1. / 30 - f * (5. / 66 - f * 691. / 2730)))))


@njit(cache=True, error_model="numpy")
def betaln(a, b):
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


@njit(cache=True, error_model="numpy")
def expit(x):
    if x >= 0:
        return 1. / (1. + math.exp(-x))
    z = math.exp(x)
    return z / (1. + z)


@njit(cache=True, error_model="numpy")
def row_terms(kind, start, stop, log_p, log_q, base, first, last, terms, row, mean, kappa):
    """log E_Beta[L_row] at mean logit m: value, d/dm, d2/dm2, kappa d/dkappa."""
    s, r = expit(mean), expit(-mean)
    a, c = max(kappa * s, 1e-300), max(kappa * r, 1e-300)
    trigamma_total = trigamma(a + c)
    if kind[row] == 1:
        x, y, constant = terms[row, 0], terms[row, 1], terms[row, 2]
        value = betaln(a + x, c + y) - betaln(a, c) + constant
        total = digamma(a + c + x + y)
        da = digamma(a + x) - total - digamma(a) + digamma(a + c)
        dc = digamma(c + y) - total - digamma(c) + digamma(a + c)
        both = -trigamma(a + c + x + y) + trigamma_total
        daa = trigamma(a + x) - trigamma(a) + both
        dcc = trigamma(c + y) - trigamma(c) + both
        dac = both
    else:
        left = first[row] - LIMIT * a - math.log(a)
        right = last[row] - LIMIT * c - math.log(c)
        peak = max(left, right)
        for i in range(start[row], stop[row]):
            peak = max(peak, base[i] + a * log_p[i] + c * log_q[i])
        mass = math.exp(left - peak) + math.exp(right - peak)
        mean_u = mean_v = square_u = square_v = cross = 0.
        for i in range(start[row], stop[row]):
            weight = math.exp(base[i] + a * log_p[i] + c * log_q[i] - peak)
            mass += weight
            mean_u += weight * log_p[i]
            mean_v += weight * log_q[i]
            square_u += weight * log_p[i] * log_p[i]
            square_v += weight * log_q[i] * log_q[i]
            cross += weight * log_p[i] * log_q[i]
        mean_u, mean_v, square_u, square_v, cross = mean_u / mass, mean_v / mass, square_u / mass, square_v / mass, cross / mass
        left_weight, right_weight = math.exp(left - peak) / mass, math.exp(right - peak) / mass
        tail_u, tail_v = -LIMIT - 1. / a, -LIMIT - 1. / c
        extra_u = extra_v = 0.
        if left_weight > 0:
            mean_u += left_weight * tail_u
            square_u += left_weight * tail_u * tail_u
            extra_u = left_weight / (a * a)
        if right_weight > 0:
            mean_v += right_weight * tail_v
            square_v += right_weight * tail_v * tail_v
            extra_v = right_weight / (c * c)
        normalizer = digamma(a + c)
        value = peak + math.log(mass) - betaln(a, c)
        da = mean_u - digamma(a) + normalizer
        dc = mean_v - digamma(c) + normalizer
        daa = square_u - mean_u * mean_u + extra_u - trigamma(a) + trigamma_total
        dcc = square_v - mean_v * mean_v + extra_v - trigamma(c) + trigamma_total
        dac = cross - mean_u * mean_v + trigamma_total
    scale = kappa * s * r
    slope = scale * (da - dc)
    curvature = scale * (r - s) * (da - dc) + scale * scale * (daa - 2 * dac + dcc)
    return value, slope, curvature, a * da + c * dc


@njit(cache=True, error_model="numpy")
def subject_terms(rows, eta, kappa, sigma, b, kind, start, stop, log_p, log_q, base, first, last, terms, d_mean, d_kappa):
    """f_u(b), d/db, d2/db2 at points b; row scores written to d_mean, d_kappa (rows x points)."""
    q = len(b)
    value = np.empty(q)
    slope = np.empty(q)
    curvature = np.empty(q)
    for j in range(q):
        value[j] = -b[j] * b[j] / (2 * sigma * sigma) - math.log(sigma) - .5 * math.log(2 * math.pi)
        slope[j] = -b[j] / (sigma * sigma)
        curvature[j] = -1. / (sigma * sigma)
    for k in range(len(rows)):
        row = rows[k]
        for j in range(q):
            integral, dm, d2m, dk = row_terms(kind, start, stop, log_p, log_q, base, first, last, terms, row, eta[row] + b[j], kappa)
            value[j] += integral
            slope[j] += dm
            curvature[j] += d2m
            d_mean[k, j] = dm
            d_kappa[k, j] = dk
    return value, slope, curvature


@njit(cache=True, error_model="numpy")
def objective(eta, kappa, sigma, group_start, group_rows, modes, roots, log_weights, kind, start, stop, log_p, log_q, base, first, last, terms):
    """Marginal log-likelihood, per-row mean scores, log-kappa and log-sigma scores.

    modes (one per subject) are read as warm starts and overwritten with the
    new posterior modes, as in the numpy GridModel.
    """
    total = 0.
    row_score = np.zeros(len(eta))
    kappa_score = 0.
    sigma_score = 0.
    one = np.empty(1)
    for u in range(len(group_start) - 1):
        rows = group_rows[group_start[u]:group_start[u + 1]]
        scratch_mean = np.empty((len(rows), 1))
        scratch_kappa = np.empty((len(rows), 1))
        b = modes[u]
        one[0] = b
        value, slope, curvature = subject_terms(rows, eta, kappa, sigma, one, kind, start, stop, log_p, log_q, base, first, last, terms, scratch_mean, scratch_kappa)
        current_value, current_slope, current_curvature = value[0], slope[0], curvature[0]
        for _ in range(100):
            step = -current_slope / current_curvature if current_curvature < 0 else current_slope * sigma * sigma
            for _ in range(40):
                one[0] = b + step
                value, slope, curvature = subject_terms(rows, eta, kappa, sigma, one, kind, start, stop, log_p, log_q, base, first, last, terms, scratch_mean, scratch_kappa)
                if value[0] >= current_value - 1e-12:
                    break
                step /= 2
            b = one[0]
            current_value, current_slope, current_curvature = value[0], slope[0], curvature[0]
            if abs(step) < 1e-9 or abs(current_slope) < 1e-9:
                break
        modes[u] = b
        scale = 1. / math.sqrt(-current_curvature) if current_curvature < 0 else sigma
        points = b + math.sqrt(2.) * scale * roots
        d_mean = np.empty((len(rows), len(roots)))
        d_kappa = np.empty((len(rows), len(roots)))
        value, _, _ = subject_terms(rows, eta, kappa, sigma, points, kind, start, stop, log_p, log_q, base, first, last, terms, d_mean, d_kappa)
        joint = log_weights + roots * roots + value
        peak = joint.max()
        evidence = peak + math.log(np.exp(joint - peak).sum())
        total += evidence + math.log(math.sqrt(2.) * scale)
        posterior = np.exp(joint - evidence)
        for k in range(len(rows)):
            row_score[rows[k]] += (posterior * d_mean[k]).sum()
            kappa_score += (posterior * d_kappa[k]).sum()
        sigma_score += (posterior * (points * points / (sigma * sigma) - 1.)).sum()
    return total, row_score, kappa_score, sigma_score


@njit(cache=True, error_model="numpy")
def profiled_precursor_log_likelihood(t, counts, precursor, target, rest, total, precursor_total, target_total, rest_total):
    """Primer log-likelihood at logits t, maximized over the precursor scale.

    Class k has probability proportional to x a_k + psi b_k + (1 - psi) d_k,
    x >= 0. With w = x A / (x A + M), A and M the precursor and mature
    opportunity totals, the class probability is w p_k + (1 - w) q_k (p, q the
    normalized precursor and mature class distributions), so the
    log-likelihood is concave in w and its maximizer is the root of a
    decreasing score, found by safeguarded Newton and bisection. w is kept in
    the range of log x in [-25, 15]. Only observed classes are passed
    (counts > 0); total and the *_total arguments sum over all classes.
    """
    result = np.empty(len(t))
    p = precursor / precursor_total
    for j in range(len(t)):
        psi = expit(t[j])
        mature_total = psi * target_total + (1 - psi) * rest_total
        q = (psi * target + (1 - psi) * rest) / mature_total
        low = math.exp(-25.) * precursor_total / (math.exp(-25.) * precursor_total + mature_total)
        high = math.exp(15.) * precursor_total / (math.exp(15.) * precursor_total + mature_total)
        w = low
        score = 0.
        for k in range(len(counts)):
            score += counts[k] * (p[k] - q[k]) / (low * p[k] + (1 - low) * q[k])
        if score > 0:
            w = high
            score = 0.
            for k in range(len(counts)):
                score += counts[k] * (p[k] - q[k]) / (high * p[k] + (1 - high) * q[k])
            if score < 0:
                w = .5 * (low + high)
                for _ in range(200):
                    score, slope = 0., 0.
                    for k in range(len(counts)):
                        mix = w * p[k] + (1 - w) * q[k]
                        difference = p[k] - q[k]
                        score += counts[k] * difference / mix
                        slope -= counts[k] * difference * difference / (mix * mix)
                    if score > 0:
                        low = w
                    else:
                        high = w
                    step = -score / slope if slope < 0 else 0.
                    proposal = w + step
                    if not (low < proposal < high):
                        proposal = .5 * (low + high)
                    if abs(proposal - w) < 1e-15 or high - low < 1e-15:
                        w = proposal
                        break
                    w = proposal
        value = 0.
        for k in range(len(counts)):
            value += counts[k] * math.log(max(w * p[k] + (1 - w) * q[k], 1e-300))
        result[j] = value
    return result
