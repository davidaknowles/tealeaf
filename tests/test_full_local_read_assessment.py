import json

import numpy as np
import pandas as pd
import pytest

from extra_scripts.assess_full_local_read_models import expected_requests, validate_shard
from extra_scripts.full_local_read_models import FIT_SETTINGS


def fixture():
    family = pd.DataFrame(dict(feature_id=['e', 'unavailable'], gene_id='g', event_type='SE'))
    contexts = pd.DataFrame([dict(gene_id='g', level_a='A', level_b='B', subjects=json.dumps(['s0', 's1', 's2', 's3']))])
    expected = expected_requests(family, contexts)
    table = expected.drop(columns='subjects').assign(converged=[True, False], p_value=[.01, 1.], raw_p_value=[.01, 1.], F_reference_p_value=[.03, 1.], statistic=[8., 0.], model_version=FIT_SETTINGS['model_version'], cohort='full', marker_variant='all', n_subjects=[4, 0], effect_size=[.5, np.nan])
    recipe = dict(shard_count=2, settings=FIT_SETTINGS, code_hashes=dict(source='frozen'))
    receipt = dict(complete=True, cohort='full', marker_variant='all', shard_index=0, shard_count=2, settings=FIT_SETTINGS, code_hashes=dict(source='frozen'), requested_tests=2, completed_tests=2, usable=1)
    return table, expected, recipe, receipt


def test_full_assessment_retains_original_unavailable_identity_and_subjects():
    table, expected, recipe, receipt = fixture()
    result = validate_shard(table, expected, recipe, receipt, 'full', 'all', 0)
    assert len(result) == 2
    assert result.loc[~result.converged, 'p_value'].eq(1.).all()


@pytest.mark.parametrize('problem', ['missing', 'duplicate', 'failed_p', 'failed_statistic', 'changed_labels', 'changed_subjects', 'partial_receipt', 'wrong_code', 'wrong_cohort', 'nan_p', 'too_few_subjects', 'nan_statistic', 'negative_statistic', 'wrong_raw_p', 'too_many_subjects'])
def test_full_assessment_refuses_partial_or_changed_hypothesis_families(problem):
    table, expected, recipe, receipt = fixture()
    if problem == 'missing':
        table = table.iloc[:1]
    elif problem == 'duplicate':
        table = pd.concat([table, table.iloc[:1]])
    elif problem == 'failed_p':
        table.loc[1, 'F_reference_p_value'] = .01
    elif problem == 'failed_statistic':
        table.loc[1, 'statistic'] = 1.
    elif problem == 'changed_labels':
        table.loc[0, 'level_a'] = 'C'
    elif problem == 'changed_subjects':
        table.loc[0, 'n_requested_subjects'] = 3
    elif problem == 'partial_receipt':
        receipt['complete'] = False
    elif problem == 'wrong_code':
        receipt['code_hashes'] = dict(source='changed')
    elif problem == 'wrong_cohort':
        table.loc[0, 'cohort'] = 'fold0'
    elif problem == 'nan_p':
        table.loc[0, 'p_value'] = np.nan
    elif problem == 'too_few_subjects':
        table.loc[0, 'n_subjects'] = 3
    elif problem == 'nan_statistic':
        table.loc[0, 'statistic'] = np.nan
    elif problem == 'negative_statistic':
        table.loc[0, 'statistic'] = -1.
    elif problem == 'wrong_raw_p':
        table.loc[0, 'raw_p_value'] = .5
    else:
        table.loc[0, 'n_subjects'] = 5
    with pytest.raises(ValueError):
        validate_shard(table, expected, recipe, receipt, 'full', 'all', 0)


def test_declared_requests_do_not_depend_on_existing_fit_rows():
    _, expected, _, _ = fixture()
    assert set(expected.test_id) == {'e|cell_type|A|B', 'unavailable|cell_type|A|B'}
    assert expected.n_requested_subjects.eq(4).all()


def test_zero_request_shard_is_valid_without_creating_hypotheses():
    table, expected, recipe, receipt = fixture()
    family = expected[['feature_id', 'gene_id', 'event_type']]
    contexts = pd.DataFrame(columns=['gene_id', 'level_a', 'level_b', 'subjects'])
    expected = expected_requests(family, contexts)
    receipt.update(requested_tests=0, completed_tests=0, usable=0)
    # Empty streamed TSV columns are object dtype when read back.
    result = validate_shard(pd.DataFrame(columns=table.columns), expected, recipe, receipt, 'full', 'all', 0)
    assert len(result) == 0
