import json
import hashlib

import numpy as np
import pandas as pd
import pytest

from extra_scripts.audit_real_pattern_read_null import DRAWS, FIT_SETTINGS, LAWS, PATTERN_NULL_VERSION, choose_parents, validate_trials


def fixture():
    family = pd.DataFrame([dict(gene_id=f'g{i}', feature_id=f'e{i}_{j}', event_type='SE') for i in range(40) for j in range(2)])
    contexts = pd.DataFrame([dict(gene_id=f'g{i}', cohort=cohort, shard_index=i % 4, level_a='A', level_b=level, median_gene_umis=100 + i, subjects=json.dumps([f's{u}' for u in range(4 + i % 20)])) for i in range(40) for cohort in ('fold0', 'fold1', 'full') for level in ('B', 'C')])
    return family, contexts


def test_parent_selection_is_gene_cohort_unique_and_independent_of_input_order():
    family, contexts = fixture()
    first = choose_parents(family, contexts, 'fold0', 2)
    second = choose_parents(family.sample(frac=1, random_state=2), contexts.sample(frac=1, random_state=3), 'fold0', 2)
    pd.testing.assert_frame_equal(first, second)
    assert not first.gene_id.duplicated().any() and first.cohort.eq('fold0').all()
    assert first.groupby(['coverage_quartile', 'subject_stratum']).size().le(2).all()
    assert first.parent_id.str.startswith('fold0|').all()
    assert first.n_requested_subjects.ge(4).all()


def trials():
    family, contexts = fixture()
    parents = choose_parents(family, contexts, 'full', 1).iloc[:1].copy()
    original = np.zeros((int(parents.n_requested_subjects.iloc[0]), 2, 2, 2), dtype=np.int64)
    digest = hashlib.sha256(original.tobytes()).hexdigest()
    parents['original_counts_sha256'] = digest
    parents['counts_json'] = json.dumps(original.tolist())
    parent = parents.iloc[0].to_dict()
    rows = [dict(parent, law=law, draw=draw, counts_sha256=digest, converged=False, p_value=1., F_reference_p_value=1., statistic=0., model_version=FIT_SETTINGS['model_version'], generator_version=PATTERN_NULL_VERSION) for law in LAWS for draw in range(DRAWS)]
    return pd.DataFrame(rows), parents


def test_null_audit_keeps_markerless_and_failed_trials():
    table, parents = trials()
    checked = validate_trials(table, parents)
    assert len(checked) == 8 and not checked.converged.any()


@pytest.mark.parametrize('problem', ['missing', 'duplicated', 'failed_p', 'nan_p', 'nonboolean', 'wrong_parent', 'wrong_cohort', 'changed_counts', 'wrong_model', 'failed_statistic', 'wrong_draw'])
def test_null_audit_refuses_changed_or_success_only_families(problem):
    table, parents = trials()
    if problem == 'missing':
        table = table.iloc[1:]
    elif problem == 'duplicated':
        table = pd.concat([table, table.iloc[:1]])
    elif problem == 'failed_p':
        table.loc[0, 'p_value'] = .001
    elif problem == 'nan_p':
        table.loc[0, 'p_value'] = np.nan
    elif problem == 'nonboolean':
        table['converged'] = 'maybe'
    elif problem == 'wrong_parent':
        table.loc[0, 'parent_id'] = 'unknown'
    elif problem == 'wrong_cohort':
        table.loc[0, 'cohort'] = 'fold0'
    elif problem == 'changed_counts':
        table.loc[0, 'original_counts_sha256'] = 'changed'
    elif problem == 'wrong_model':
        table.loc[0, 'model_version'] = 'unknown'
    elif problem == 'failed_statistic':
        table.loc[0, 'statistic'] = 2.
    else:
        table.loc[0, 'counts_sha256'] = 'different simulation'
    with pytest.raises(ValueError):
        validate_trials(table, parents)
