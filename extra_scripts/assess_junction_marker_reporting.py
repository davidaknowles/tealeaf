#!/usr/bin/env python3
"""All-marker tests with junction-only directions on complete fixed families."""

import argparse
import json
from pathlib import Path

import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.full_local_read_models import COHORTS
from extra_scripts.assess_full_local_read_models import load_cohort, split_table
from extra_scripts.assess_paired_inference_audit import split_assessment
from tealeaf.sc.replication_audit import reexpress_event_directions, rank_direction_table, ranked_direction_summary, ranked_category_summary


def aligned_reporting_table(tests, reports, column):
    """Align complete event requests, not just successful junction fits.

    This is a dataset-driver provenance check, not a new statistical method.
    The original tests/p-values are unchanged and unavailable effects stay NaN.
    """
    keys = ['test_id', 'feature_id', 'gene_id', 'event_type', 'level_a', 'level_b', 'contrast_id', 'cohort', 'n_requested_subjects']
    if tests.test_id.duplicated().any() or reports.test_id.duplicated().any() or set(tests.test_id) != set(reports.test_id):
        raise ValueError('identical complete requested families required for reporting')
    check = tests[keys].merge(reports[keys], on='test_id', validate='one_to_one', suffixes=('_test', '_report'))
    if any(not check[f'{key}_test'].eq(check[f'{key}_report']).all() for key in keys[1:]):
        raise ValueError('original event, contrast, cohort and requested subjects must match')
    result = tests.copy()
    result['junction_reporting_effect'] = result.test_id.map(reports.set_index('test_id')[column])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--all-marker-assessment', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError('new complete reporting assessment required')
    recipe_path = args.root / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    baseline_manifest_path = args.all_marker_assessment / 'manifest.json'
    baseline_manifest = json.loads(baseline_manifest_path.read_text())
    if baseline_manifest['recipe_sha256'] != file_hash(recipe_path) or baseline_manifest['recipe'] != recipe:
        raise ValueError('completed all-marker endpoint from the original recipe required')
    tables, receipts = {}, {}
    for variant in ('all', 'junction'):
        tables[variant], receipts[variant] = {}, {}
        for cohort in COHORTS:
            tables[variant][cohort], receipts[variant][cohort] = load_cohort(args.root, recipe, cohort, variant)
    if baseline_manifest['cohorts'] != receipts['all']:
        raise ValueError('baseline LR mapping must come from these complete all-marker fits')
    mapping_path = args.all_marker_assessment / 'lr_mapping.tsv.gz'
    mapped = pd.read_csv(mapping_path, sep='\t')
    valid = mapped.mapping_complete.astype(str).str.lower().eq('true') & mapped.minimum_pooled_depth.ge(20) & mapped.pooled_replicated.notna()
    eligible = mapped.loc[valid].copy()
    eligible['pooled_replicated'] = eligible.pooled_replicated.astype(str).str.lower().eq('true')
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True)
    summaries = []
    for reporting, column in (('junction_log_odds', 'effect_size'), ('junction_pooled_fraction', 'pooled_marker_effect')):
        replaced = {cohort: aligned_reporting_table(tables['all'][cohort], tables['junction'][cohort], column) for cohort in COHORTS}
        for tail, p_column in (('native_chi1', 'p_value'), ('F_reference', 'F_reference_p_value')):
            output = args.output_dir / reporting / tail
            output.mkdir(parents=True)
            method = f'All-marker local read mixed tests, {tail}, {reporting} reporting only'
            folds = [split_table(replaced[cohort], method, p_column, 'junction_reporting_effect') for cohort in ('fold0', 'fold1')]
            split_assessment(folds, repo, output, method)
            baseline = eligible.copy()
            if tail == 'F_reference':
                pvalues = tables['all']['full'].set_index(['feature_id', 'contrast_id'])[p_column]
                baseline['p_value'] = [pvalues.loc[(row.feature_id, row.contrast_id)] for row in baseline.itertuples(index=False)]
                baseline['raw_p_value'] = baseline.p_value
            reported = reexpress_event_directions(baseline, replaced['full'], 'junction_reporting_effect', method)
            original_rank = rank_direction_table(baseline, len(baseline))
            ranked = rank_direction_table(reported, len(reported))
            keys = ['feature_id', 'contrast_id', 'rank']
            if original_rank[keys].to_records(index=False).tolist() != ranked[keys].to_records(index=False).tolist():
                raise ValueError('junction reporting cannot change eligible LR identities or test ranks')
            summary = pd.DataFrame(ranked_direction_summary(ranked))
            summary.to_csv(output / 'lr_rank_summary.tsv', sep='\t', index=False)
            pd.DataFrame(ranked_category_summary(ranked)).to_csv(output / 'lr_event_type_composition.tsv', sep='\t', index=False)
            ranked.head(200).to_csv(output / 'lr_rank.tsv.gz', sep='\t', index=False)
            summaries.append(summary.assign(reporting=reporting, tail=tail))
    pd.concat(summaries, ignore_index=True).to_csv(args.output_dir / 'lr_summary.tsv', sep='\t', index=False)
    manifest = dict(recipe_sha256=file_hash(recipe_path), recipe=recipe, receipts=receipts, input_hashes={str(path): file_hash(path) for path in (baseline_manifest_path, mapping_path)}, code_hashes={str(path): file_hash(path) for path in (Path(__file__), Path(reexpress_event_directions.__code__.co_filename), Path(split_assessment.__code__.co_filename))}, scope='junction-only reporting for complete all-marker tests, unchanged hypotheses, p-values, significance ranks, split universes and LR-evaluable family; missing/zero LR replacement directions count as nonagreements', limitation='marker reporting is not absolute RNA PSI or correction for outside event classes/type-specific capture; reporting cannot repair test calibration, split effect agreement must also be assessed', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
