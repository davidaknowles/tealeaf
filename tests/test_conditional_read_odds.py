import numpy as np
import pytest
from scipy import integrate, special, stats

from tealeaf.sc.conditional_read_odds import ConditionalReadOdds, conditional_read_odds_test, simulate_local_read_counts, simulate_conditional_read_counts


def panel():
    return np.array([[[[8, 12], [11, 9]], [[5, 4], [3, 8]]], [[[2, 10], [6, 4]], [[1, 7], [2, 7]]], [[[5, 5], [7, 3]], [[4, 5], [8, 2]]], [[[10, 3], [8, 6]], [[7, 8], [4, 9]]]])


def test_exact_conditional_stratum_matches_noncentral_hypergeometric():
    likelihood = ConditionalReadOdds(panel())
    table = panel()[0, 0]
    first, second = table.sum(axis=1)
    included = table[:, 0].sum()
    effect = np.array([-3., -.2, 0., .7, 3.])
    values, score, info = likelihood.subject_strata[0][0].evaluate(effect)
    expected = stats.nchypergeom_fisher.logpmf(table[1, 0], first + second, included, second, np.exp(effect))
    np.testing.assert_allclose(values, expected, atol=1e-12)
    step = 1e-5
    upper = likelihood.subject_strata[0][0].evaluate(effect + step)
    lower = likelihood.subject_strata[0][0].evaluate(effect - step)
    np.testing.assert_allclose((upper[0] - lower[0]) / (2 * step), score, atol=1e-8)
    np.testing.assert_allclose(-(upper[1] - lower[1]) / (2 * step), info, atol=1e-8)


@pytest.mark.parametrize("sd", [0., .01, .6, 2., 5.])
def test_shared_subject_integral_matches_adaptive_numeric_integration(sd):
    likelihood = ConditionalReadOdds(panel())
    mean = .7
    expected = 0.
    for index in range(4):
        if sd == 0:
            expected -= likelihood.subject_evaluate(index, [mean])[0][0]
        else:
            value = integrate.quad(lambda z: np.exp(likelihood.subject_evaluate(index, [mean + sd * z])[0][0]) * stats.norm.pdf(z), -12., 12., epsabs=1e-12, epsrel=1e-10)[0]
            expected -= np.log(value)
    actual = likelihood.marginal_objective(mean, sd, nodes=61, gradient=False)
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-8)


def test_adaptive_quadrature_derivatives_match_objective_differences():
    likelihood = ConditionalReadOdds(panel())
    mean, sd, step = .3, .8, 1e-5
    _, derivative = likelihood.marginal_objective(mean, sd, nodes=41)
    values = [(likelihood.marginal_objective(mean + step, sd, nodes=41, gradient=False) - likelihood.marginal_objective(mean - step, sd, nodes=41, gradient=False)) / (2 * step), (likelihood.marginal_objective(mean, sd * np.exp(step), nodes=41, gradient=False) - likelihood.marginal_objective(mean, sd * np.exp(-step), nodes=41, gradient=False)) / (2 * step)]
    np.testing.assert_allclose(derivative, values, atol=2e-7)


def test_type_swap_reverses_effect_and_preserves_marginal_likelihood():
    original = ConditionalReadOdds(panel())
    swapped = ConditionalReadOdds(panel()[:, :, ::-1, :])
    for mean in (-3., 0., 2.):
        first, first_score = original.marginal_objective(mean, .7, nodes=41)
        second, second_score = swapped.marginal_objective(-mean, .7, nodes=41)
        np.testing.assert_allclose(first, second, atol=1e-10)
        np.testing.assert_allclose(first_score, second_score * [-1, 1], atol=1e-10)


def test_primer_order_and_uninformative_strata_do_not_change_likelihood():
    counts = panel()
    original = ConditionalReadOdds(counts)
    reversed_primers = ConditionalReadOdds(counts[:, ::-1])
    added = np.zeros((4, 1, 2, 2), dtype=int)
    added[:, 0, :, 0] = [3, 50]
    extra = ConditionalReadOdds(np.concatenate([counts, added], axis=1))
    value = original.marginal_objective(.4, .7, nodes=41)
    for other in (reversed_primers, extra):
        actual = other.marginal_objective(.4, .7, nodes=41)
        np.testing.assert_allclose(actual[0], value[0], atol=1e-12)
        np.testing.assert_allclose(actual[1], value[1], atol=1e-12)
    np.testing.assert_array_equal(counts, panel())


@pytest.mark.parametrize("counts", [np.ones((2, 2)), np.full((4, 2, 2, 2), .5), np.full((4, 2, 2, 2), -1), np.full((4, 2, 2, 2), np.nan)])
def test_invalid_counts_fail(counts):
    with pytest.raises(ValueError):
        ConditionalReadOdds(counts)


def test_requested_unmeasured_subject_does_not_create_information():
    counts = np.concatenate([np.zeros((1, 2, 2, 2), dtype=int), panel()])
    likelihood = ConditionalReadOdds(counts)
    assert likelihood.informative_subjects.tolist() == [1, 2, 3, 4]
    assert likelihood.n_informative_subjects == 4
    with pytest.raises(ValueError, match="four"):
        conditional_read_odds_test(ConditionalReadOdds(counts[:4]))


def test_nested_fit_and_type_reversal_preserve_test_and_effect_direction():
    counts = np.repeat(np.array([[[[15, 35], [30, 20]], [[8, 22], [18, 12]]]]), 8, axis=0)
    fitted = conditional_read_odds_test(ConditionalReadOdds(counts), nodes=21)
    swapped = conditional_read_odds_test(ConditionalReadOdds(counts[:, :, ::-1, :]), nodes=21)
    assert fitted["converged"] and swapped["converged"]
    assert fitted["p_value"] < 1e-6
    assert fitted["log_odds_effect"] > 0
    assert swapped["log_odds_effect"] < 0
    np.testing.assert_allclose(fitted["statistic"], swapped["statistic"], rtol=1e-7)
    assert fitted["alternative_objective"] <= fitted["null_objective"]
    assert fitted["quadrature_error"] < 1e-4


def test_extreme_boundary_observation_has_finite_conditional_objective():
    counts = np.repeat(np.array([[[[0, 2000], [1, 999]]]]), 4, axis=0)
    likelihood = ConditionalReadOdds(counts)
    for mean, sd in ((0., 0.), (20., 3.), (-20., 3.)):
        value, derivative = likelihood.marginal_objective(mean, sd, nodes=41)
        assert np.isfinite(value) and np.isfinite(derivative).all()


def test_unconditional_generator_preserves_type_totals_and_random_stream():
    totals = panel().sum(axis=-1)
    baselines = np.array([[-1., 1.]] * 4)
    effects = np.array([-3., 0., 1., 3.])
    first = simulate_local_read_counts(totals, baselines, effects, np.random.default_rng(947))
    second = simulate_local_read_counts(totals, baselines, effects, np.random.default_rng(947))
    np.testing.assert_array_equal(first.sum(axis=-1), totals)
    np.testing.assert_array_equal(first, second)
    assert np.issubdtype(first.dtype, np.integer) and (first >= 0).all()


def test_conditional_generator_preserves_both_margins_and_zero_strata():
    counts = np.concatenate([panel(), np.zeros((1, 2, 2, 2), dtype=int)])
    counts[-1, 0] = [[3, 0], [4, 0]]
    original = counts.copy()
    generated = simulate_conditional_read_counts(counts, [.4, -.7, 0., 1., 0.], np.random.default_rng(811))
    np.testing.assert_array_equal(generated.sum(axis=-1), counts.sum(axis=-1))
    np.testing.assert_array_equal(generated.sum(axis=-2), counts.sum(axis=-2))
    np.testing.assert_array_equal(generated[-1], counts[-1])
    np.testing.assert_array_equal(counts, original)


@pytest.mark.parametrize("problem", ["fractional_totals", "wrong_baseline", "wrong_effects", "nan_baseline"])
def test_binomial_generator_rejects_invalid_geometry(problem):
    totals = panel().sum(axis=-1).astype(float)
    baseline, effects = np.zeros((4, 2)), np.zeros(4)
    if problem == "fractional_totals":
        totals[0, 0, 0] += .5
    elif problem == "wrong_baseline":
        baseline = baseline[:, :1]
    elif problem == "wrong_effects":
        effects = effects[:3]
    else:
        baseline[0, 0] = np.nan
    with pytest.raises(ValueError):
        simulate_local_read_counts(totals, baseline, effects, np.random.default_rng(11))
