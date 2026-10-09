import json
import sys

import pandas as pd
import pytest

from extra_scripts.collate_full_local_read_counts import summarize_markers, gene_key_counts
from extra_scripts import collate_full_local_read_counts as driver
from extra_scripts.audit_event_local_read_support import file_hash


def signatures():
    return pd.DataFrame(dict(feature_id=['e'] * 6, subject='s', cell_type='A', primer='RH', signature=[0, 17, 8, 25, 5, 10], count=[11, 2, 3, 7, 4, 5]))


def test_local_and_junction_ablation_count_keys_once_and_exclude_conflicts():
    table = summarize_markers(signatures())
    assert table.iloc[0][['included', 'excluded', 'junction_included', 'junction_excluded', 'gene_keys']].tolist() == [6, 8, 4, 15, 32]


def test_signature_aggregation_requires_exact_unique_keys():
    source = signatures()
    with pytest.raises(ValueError, match='unique'):
        summarize_markers(pd.concat([source, source.iloc[:1]]))
    for column, value in (('count', -.1), ('count', .5), ('signature', 32)):
        invalid = source.astype({column: float}).copy()
        invalid.loc[0, column] = value
        with pytest.raises(ValueError, match='exact'):
            summarize_markers(invalid)


def test_gene_molecule_totals_are_not_summed_across_events_and_zero_genes_remain_declared():
    marker = summarize_markers(signatures())
    marker = pd.concat([marker, marker.assign(feature_id='f')])
    family = pd.DataFrame(dict(feature_id=['e', 'f', 'no_reads'], gene_id=['g', 'g', 'zero']))
    genes = gene_key_counts(marker, family)
    assert len(genes) == 1
    assert genes.iloc[0].gene_keys == 32
    assert len(family) == 3
    with pytest.raises(ValueError, match='identical'):
        gene_key_counts(marker.iloc[:1], family)
    invalid = marker.copy()
    invalid.iloc[0, invalid.columns.get_loc('gene_keys')] = 31
    with pytest.raises(ValueError, match='identical'):
        gene_key_counts(invalid, family)


def whole_fixture(tmp_path, monkeypatch, problem=None):
    source, output = tmp_path / 'source', tmp_path / 'output'
    source.mkdir()
    groups = {f'barcode{i}': ['s', 'A', 'RH'] for i in range(8)}
    packets = [dict(library=f'lib{i}', barcode_groups={f'barcode{i}': groups[f'barcode{i}']}) for i in range(8)]
    recipe = dict(full_catalog=True, library_union=True, production_cell_qc=dict(exact_cached_group_and_primer_total_match=True), bams=packets, barcode_groups=groups, declared_events=['e', 'f', 'unavailable'], events=[dict(event_id=key) for key in ('e', 'f')], selection='whole annotation family')
    recipe_path = source / 'recipe.json'
    recipe_path.write_text(json.dumps(recipe))
    family = pd.DataFrame(dict(feature_id=['e', 'f', 'unavailable'], gene_id=['g', 'g', 'zero'], event_type='SE'))
    family.to_csv(source / 'declared_family.tsv.gz', sep='\t', index=False)
    family.assign(status=['available', 'available', 'unavailable']).to_csv(source / 'feature_definitions.tsv', sep='\t', index=False)
    for index, packet in enumerate(packets):
        folder = source / f'shard_{index}'
        folder.mkdir()
        counts = pd.concat([signatures(), signatures().assign(feature_id='f')])
        counts.to_csv(folder / 'support.tsv.gz', sep='\t', index=False)
        receipt = dict(input=packet, complete=True, shard=index, recipe_sha256=file_hash(recipe_path), requested_events=2, filters=dict(barcode_umi_event_keys=int(counts['count'].sum())))
        if index == 7 and problem == 'partial':
            receipt['complete'] = False
        elif index == 7 and problem == 'wrong_total':
            receipt['filters']['barcode_umi_event_keys'] -= 1
        elif index == 7 and problem == 'changed_recipe':
            receipt['recipe_sha256'] = 'changed'
        (folder / 'manifest.json').write_text(json.dumps(receipt))
    monkeypatch.setattr(sys, 'argv', ['collate', '--source', str(source), '--output-dir', str(output), '--shard-count', '2'])
    return output


def test_whole_count_collator_keeps_markerless_definitions_and_all_libraries(tmp_path, monkeypatch):
    output = whole_fixture(tmp_path, monkeypatch)
    driver.main()
    manifest = json.loads((output / 'manifest.json').read_text())
    assert manifest['declared_events'] == 3
    assert manifest['available_marker_geometry'] == 2
    assert len(manifest['library_receipts']) == 8
    assert manifest['production_barcodes'] == 8
    family = pd.read_csv(output / 'declared_family.tsv.gz', sep='\t')
    assert set(family.feature_id) == {'e', 'f', 'unavailable'}
    genes = pd.read_csv(output / 'gene_keys.tsv.gz', sep='\t')
    assert len(genes) == 1 and genes.iloc[0].gene_keys == 256
    counts = pd.read_csv(output / 'shard_0/counts.tsv.gz', sep='\t')
    assert len(counts) == 2 and counts.included.eq(48).all()


@pytest.mark.parametrize('problem', ['partial', 'wrong_total', 'changed_recipe'])
def test_whole_count_collator_rejects_incompatible_or_incomplete_source(tmp_path, monkeypatch, problem):
    output = whole_fixture(tmp_path, monkeypatch, problem)
    with pytest.raises(ValueError):
        driver.main()
    assert not output.exists()
