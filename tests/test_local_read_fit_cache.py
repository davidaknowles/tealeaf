import numpy as np
import pytest

from tealeaf.sc.local_read_fit_cache import ExactLocalReadFitCache


def counts():
    return np.arange(32, dtype=np.int64).reshape(4, 2, 2, 2)


def test_exact_cache_reuses_counts_not_outcomes_and_preserves_result_values():
    calls = []
    def fit(values):
        calls.append(values.copy())
        return dict(p_value=.02, effect=-.5, nested=[1, 2])
    cache = ExactLocalReadFitCache(fit)
    first, error, hit = cache.evaluate(counts())
    assert not hit and not error
    first['nested'].append(9)
    second, error, hit = cache.evaluate(counts().copy())
    assert hit and not error and second == dict(p_value=.02, effect=-.5, nested=[1, 2])
    assert len(calls) == cache.evaluations == cache.hits == 1


def test_same_bytes_different_shape_or_reordered_counts_cannot_share_fit():
    cache = ExactLocalReadFitCache(lambda value: dict(shape=value.shape, first=int(value.ravel()[0])))
    for value in (counts(), counts().reshape(2, 4, 2, 2), counts()[::-1]):
        _, _, hit = cache.evaluate(value)
        assert not hit
    assert cache.evaluations == 3


def test_failures_remain_failures_and_only_numerical_exceptions_are_cached():
    def failure(values):
        raise ValueError('insufficient evidence')
    cache = ExactLocalReadFitCache(failure)
    assert cache.evaluate(counts()) == (None, 'insufficient evidence', False)
    assert cache.evaluate(counts()) == (None, 'insufficient evidence', True)
    with pytest.raises(RuntimeError):
        ExactLocalReadFitCache(lambda values: (_ for _ in ()).throw(RuntimeError('external failure'))).evaluate(counts())


def test_cache_is_bounded_and_does_not_reuse_evicted_or_independent_recipes():
    cache = ExactLocalReadFitCache(lambda value: dict(p=1), max_entries=1)
    for value in (counts(), counts() + 1, counts()):
        assert cache.evaluate(value)[2] is False
    assert len(cache.entries) == 1 and cache.evaluations == 3
    assert ExactLocalReadFitCache(lambda value: dict(p=.1)).evaluate(counts())[0]['p'] == .1


def test_real_marginal_fit_is_unchanged_by_exact_reuse():
    from extra_scripts.full_local_read_models import marginal_fit
    original = np.array([[[[2, 3], [4, 2]], [[3, 4], [4, 3]]], [[[1, 4], [2, 4]], [[2, 4], [3, 5]]], [[[3, 5], [4, 4]], [[1, 6], [2, 4]]], [[[4, 3], [2, 3]], [[2, 1], [2, 4]]]], dtype=np.int64)
    uncached = marginal_fit(original)
    cache = ExactLocalReadFitCache(marginal_fit)
    first, error, hit = cache.evaluate(original)
    second, error_again, hit_again = cache.evaluate(original)
    assert not error and not error_again and not hit and hit_again
    assert first.keys() == second.keys() == uncached.keys()
    for key, value in uncached.items():
        if isinstance(value, (float, np.floating)):
            np.testing.assert_allclose([first[key], second[key]], value, rtol=0, atol=0, equal_nan=True)
        else:
            assert first[key] == second[key] == value


@pytest.mark.parametrize('value', [counts().astype(float), counts() - 1, counts()[:, :, 0], np.full((1, 1, 2, 2), 2 ** 51, dtype=np.int64)])
def test_invalid_original_counts_are_not_repaired_or_used_as_cache_keys(value):
    cache = ExactLocalReadFitCache(lambda value: dict(p=1))
    with pytest.raises(ValueError):
        cache.evaluate(value)
