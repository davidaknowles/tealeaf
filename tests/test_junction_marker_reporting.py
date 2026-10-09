import json
import sys

import numpy as np
import pandas as pd
import pytest

from extra_scripts.assess_junction_marker_reporting import aligned_reporting_table
from extra_scripts import assess_junction_marker_reporting as driver
from extra_scripts.audit_event_local_read_support import file_hash


def fixture():
    rows = [dict(test_id=f't{i}', feature_id=f'e{i}', gene_id='g', event_type='SE', level_a='A', level_b='B', contrast_id='A_B', cohort='full', n_requested_subjects=8, effect_size=.4, pooled_marker_effect=.1, p_value=.001) for i in range(2)]
    tests = pd.DataFrame(rows)
    reports = tests.iloc[::-1].copy()
    reports.loc[reports.test_id.eq('t0'), 'effect_size'] = np.nan
    reports.loc[reports.test_id.eq('t1'), 'effect_size'] = -.2
    return tests, reports


def test_reporting_preserves_test_values_order_and_failed_effects():
    tests, reports = fixture()
    result = aligned_reporting_table(tests, reports, 'effect_size')
    pd.testing.assert_frame_equal(result[tests.columns], tests)
    assert np.isnan(result.junction_reporting_effect.iloc[0])
    assert result.junction_reporting_effect.iloc[1] == -.2
    assert len(result) == len(tests)


@pytest.mark.parametrize('problem', ['missing', 'duplicate', 'subject', 'type', 'cohort'])
def test_partial_or_changed_reporting_families_cannot_be_aligned(problem):
    tests, reports = fixture()
    if problem == 'missing':
        reports = reports.iloc[:1]
    elif problem == 'duplicate':
        reports = pd.concat([reports, reports.iloc[:1]])
    elif problem == 'subject':
        reports.loc[reports.index[0], 'n_requested_subjects'] = 7
    elif problem == 'type':
        reports.loc[reports.index[0], 'level_a'] = 'C'
    else:
        reports.loc[reports.index[0], 'cohort'] = 'fold0'
    with pytest.raises(ValueError):
        aligned_reporting_table(tests, reports, 'effect_size')


def test_reporting_driver_keeps_all_marker_pvalues_and_fixed_LR_ranks(tmp_path, monkeypatch):
    root, baseline, output = [tmp_path / name for name in ('root', 'baseline', 'output')]
    root.mkdir()
    baseline.mkdir()
    recipe_path = root / 'recipe.json'
    recipe_path.write_text(json.dumps(dict(original='recipe')))
    receipts = {cohort: [dict(cohort=cohort, marker_variant='all')] for cohort in driver.COHORTS}
    (baseline / 'manifest.json').write_text(json.dumps(dict(recipe_sha256=file_hash(recipe_path), recipe=dict(original='recipe'), cohorts=receipts)))
    tests, reports = fixture()
    tests['raw_p_value'] = tests.p_value
    tests['F_reference_p_value'] = [.01, .02]
    tests['statistic'] = [10., 8.]
    tests['converged'] = True
    tests['median_gene_umis'] = 100.
    reports['p_value'] = .9
    mapped = tests[['feature_id', 'gene_id', 'event_type', 'contrast_id', 'p_value', 'raw_p_value', 'statistic']].assign(method='original', mapping_complete=True, minimum_pooled_depth=30., short_read_effect=.4, long_read_effect=.2, pooled_replicated=True, replicate_1_dot_product=.08, replicate_2_dot_product=.08)
    mapped.to_csv(baseline / 'lr_mapping.tsv.gz', sep='\t', index=False)
    def load(root, recipe, cohort, variant):
        return (tests if variant == 'all' else reports).assign(cohort=cohort), [dict(cohort=cohort, marker_variant=variant)]
    calls = []
    def split(folds, repo, output, method):
        for fold in folds:
            assert fold.p_value.tolist() in ([.001, .001], [.01, .02])
            assert fold.test_id.tolist() == tests.test_id.tolist()
        calls.append(method)
    monkeypatch.setattr(driver, 'load_cohort', load)
    monkeypatch.setattr(driver, 'split_assessment', split)
    monkeypatch.setattr(sys, 'argv', ['driver', '--root', str(root), '--all-marker-assessment', str(baseline), '--output-dir', str(output)])
    driver.main()
    assert len(calls) == 4
    for tail in ('native_chi1', 'F_reference'):
        ranked = pd.read_csv(output / 'junction_log_odds' / tail / 'lr_rank.tsv.gz', sep='\t')
        assert ranked.feature_id.tolist() == ['e0', 'e1']
        assert not ranked.pooled_replicated.any()
        assert ranked.direction_available.tolist() == [False, True]
        pooled = pd.read_csv(output / 'junction_pooled_fraction' / tail / 'lr_rank.tsv.gz', sep='\t')
        assert pooled.feature_id.tolist() == ranked.feature_id.tolist()
        assert pooled.pooled_replicated.all()
