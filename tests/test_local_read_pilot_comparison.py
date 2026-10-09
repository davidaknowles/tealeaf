import pandas as pd
import pytest

from extra_scripts.compare_local_read_pilots import matched_pilots


def fixture():
    cases = pd.DataFrame([dict(fold=0, test_id='t')])
    frames = []
    for version in ('v1', 'v2'):
        frames.append(pd.DataFrame([dict(fold=0, test_id='t', model=model, variant=variant, p_value=1., converged=False, counts_sha256=variant, requested_subjects=4, model_version=f'local_read_binomial_random_intercept_slope_{version}' if model == 'unconditional' else 'conditional_local_read_random_slope_v1') for model in ('conditional', 'unconditional') for variant in ('all local markers', 'junction markers')]))
    return *frames, cases


def test_matched_pilots_keep_every_failed_identity_and_original_counts():
    old, new, cases = fixture()
    table = matched_pilots(old, new, cases)
    assert len(table) == 4 and not table.converged_new.any()


def test_same_model_adaptive_comparison_requires_explicit_versions():
    old, new, cases = fixture()
    old.loc[old.model.eq('unconditional'), 'model_version'] = 'local_read_binomial_random_intercept_slope_v2'
    with pytest.raises(ValueError):
        matched_pilots(old, new, cases)
    assert len(matched_pilots(old, new, cases, ('v2', 'v2'))) == 4


@pytest.mark.parametrize('problem', ['changed_counts', 'changed_subjects', 'changed_version'])
def test_pilot_comparison_refuses_nonmatched_model_recipes(problem):
    old, new, cases = fixture()
    if problem == 'changed_counts':
        new['counts_sha256'] = 'changed'
    elif problem == 'changed_subjects':
        new['requested_subjects'] = 3
    else:
        new.loc[new.model.eq('unconditional'), 'model_version'] = 'unrecognized'
    with pytest.raises(ValueError):
        matched_pilots(old, new, cases)
