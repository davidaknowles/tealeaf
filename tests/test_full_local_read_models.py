import argparse
import json

import numpy as np
import pandas as pd
import pytest

from extra_scripts import full_local_read_models as driver
from extra_scripts.audit_event_local_read_support import file_hash
from tealeaf.sc.local_read_support import local_read_count_tensor


def metadata_fixture():
    groups = pd.DataFrame([dict(subject=f's{i}', cell_type=level, primer=primer) for i in range(8) for level in ('A', 'B') for primer in driver.PRIMERS])
    folds = pd.DataFrame(dict(subject=[f's{i}' for i in range(8)], condition='control', fold=[0] * 4 + [1] * 4))
    return groups, folds


def full_fixture(tmp_path):
    source, root = tmp_path / 'counts', tmp_path / 'models'
    source.mkdir()
    groups, folds = metadata_fixture()
    groups.to_csv(source / 'groups.tsv', sep='\t', index=False)
    fold_path = tmp_path / 'folds.tsv'
    folds.to_csv(fold_path, sep='\t', index=False)
    family = pd.DataFrame(dict(feature_id=['SUPPA2:e', 'SUPPA2:f', 'unavailable', 'zero_gene'], gene_id=['g', 'g', 'g', 'z'], event_type='SE', shard_index=[0, 0, 0, 1]))
    family.to_csv(source / 'declared_family.tsv.gz', sep='\t', index=False)
    family.assign(status=['ok', 'ok', 'unavailable', 'ok']).to_csv(source / 'feature_definitions.tsv', sep='\t', index=False)
    genes = groups.assign(gene_id='g', gene_keys=20)
    genes.to_csv(source / 'gene_keys.tsv.gz', sep='\t', index=False)
    counts = pd.concat([groups.assign(feature_id=event, included=5, excluded=5, junction_included=2, junction_excluded=3, gene_keys=20) for event in ('SUPPA2:e', 'SUPPA2:f')])
    shards = []
    for index in range(2):
        folder = source / f'shard_{index}'
        folder.mkdir()
        local = counts if index == 0 else counts.iloc[:0]
        path = folder / 'counts.tsv.gz'
        local.to_csv(path, sep='\t', index=False)
        shards.append(dict(shard_index=index, count_rows=len(local), counts_sha256=file_hash(path)))
    manifest = dict(count_shards=shards, library_receipts=[dict(complete=True) for _ in range(8)], production_cell_qc=dict(exact_cached_group_and_primer_total_match=True), declared_events=4)
    (source / 'manifest.json').write_text(json.dumps(manifest))
    return argparse.Namespace(source=source, output_root=root, subject_folds=fold_path, cohort='full', shard_index=0, marker_variant='all')


def test_metadata_keeps_whole_cohort_or_original_fixed_subject_halves():
    groups, folds = metadata_fixture()
    assert len(driver.cohort_metadata(groups, folds, 'full')) == 16
    assert set(driver.cohort_metadata(groups, folds, 'fold0').subject) == {f's{i}' for i in range(4)}
    assert set(driver.cohort_metadata(groups, folds, 'fold1').subject) == {f's{i}' for i in range(4, 8)}
    with pytest.raises(ValueError, match='both primers'):
        driver.cohort_metadata(groups.iloc[1:], folds, 'full')


def test_coverage_contexts_require_both_types_in_original_subjects_not_marker_fits():
    groups, folds = metadata_fixture()
    metadata = driver.cohort_metadata(groups, folds, 'full')
    counts = groups.assign(gene_keys=20)
    context = driver.gene_contexts(metadata, counts)
    assert len(context) == 1 and len(json.loads(context[0]['subjects'])) == 8
    counts.loc[counts.subject.isin(['s0', 's1', 's2', 's3', 's4']) & counts.cell_type.eq('B'), 'gene_keys'] = 0
    assert driver.gene_contexts(metadata, counts) == []


def test_full_recipe_keeps_unavailable_events_and_counts_zero_gene_as_uncovered(tmp_path):
    args = full_fixture(tmp_path)
    driver.prepare(args)
    recipe = json.loads((args.output_root / 'recipe.json').read_text())
    assert recipe['declared_events'] == 4
    assert recipe['requested_tests_by_cohort'] == {'fold0': 3, 'fold1': 3, 'full': 3}
    family = pd.read_csv(args.output_root / 'declared_family.tsv.gz', sep='\t')
    assert set(family.feature_id) == {'SUPPA2:e', 'SUPPA2:f', 'unavailable', 'zero_gene'}


def test_full_worker_streams_all_requests_with_markerless_failures_at_one(tmp_path, monkeypatch):
    args = full_fixture(tmp_path)
    driver.prepare(args)
    def fit(likelihood, **kwargs):
        return dict(converged=True, p_value=.001, statistic=12., n_subjects=likelihood.n_paired_subjects, log_odds_effect=.4, model_version=driver.local_read_mixed.MODEL_VERSION)
    monkeypatch.setattr(driver.local_read_mixed, 'local_read_mixed_adaptive_test', fit)
    driver.fit(args)
    folder = args.output_root / 'all/full/shard_0'
    table = pd.read_csv(folder / 'tests.tsv.gz', sep='\t')
    assert len(table) == 3
    unavailable = table.loc[table.feature_id.eq('unavailable')].iloc[0]
    assert not unavailable.converged and unavailable.p_value == 1. and unavailable.n_requested_subjects == 8
    assert unavailable.n_local_included_keys == unavailable.n_local_excluded_keys == 0
    assert table.loc[table.converged, 'effect_size'].eq(.4).all()
    receipt = json.loads((folder / 'manifest.json').read_text())
    assert receipt['complete'] and receipt['requested_tests'] == receipt['completed_tests'] == 3
    with pytest.raises(ValueError, match='overwrite'):
        driver.fit(args)


@pytest.mark.parametrize('changed', ['family', 'contexts', 'code'])
def test_full_worker_rejects_changed_declared_scientific_recipe(tmp_path, changed):
    args = full_fixture(tmp_path)
    driver.prepare(args)
    if changed == 'code':
        path = args.output_root / 'recipe.json'
        recipe = json.loads(path.read_text())
        recipe['code_hashes'] = {'changed': 'changed'}
        path.write_text(json.dumps(recipe))
    else:
        path = args.output_root / ('declared_family.tsv.gz' if changed == 'family' else 'contexts.tsv.gz')
        path.write_bytes(b'changed')
    with pytest.raises(ValueError):
        driver.fit(args)
    assert not (args.output_root / 'all/full/shard_0').exists()


def test_generic_count_tensor_preserves_original_subject_and_primer_order():
    lookup = {('e', 's2', 'B', 'RH'): dict(included=3, excluded=7)}
    counts = local_read_count_tensor('e', ['s2', 's1'], ['A', 'B'], ['DT', 'RH'], lookup)
    assert counts.shape == (2, 2, 2, 2)
    np.testing.assert_array_equal(counts[0, 1, 1], [3, 7])
    assert not counts[1].any() and not counts[:, 0].any()
    with pytest.raises(ValueError, match='exact'):
        local_read_count_tensor('e', ['s2'], ['A', 'B'], ['RH'], {('e', 's2', 'B', 'RH'): dict(included=.2, excluded=0)})
