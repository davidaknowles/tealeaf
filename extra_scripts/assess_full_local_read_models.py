#!/usr/bin/env python3
"""Complete split and own-ranked LR endpoints for whole local-marker fits."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.full_local_read_models import COHORTS, FIT_SETTINGS, requested_id
from extra_scripts.assess_paired_inference_audit import split_assessment
from extra_scripts.evaluate_suppa2_statistics import normalize_pairs
from tealeaf.sc.replication_audit import coverage_correlation, ranked_direction_summary, ranked_category_summary, rank_direction_table, reexpress_event_directions


def expected_requests(family, contexts):
    """Original event/contrast family, never reconstructed from successful fits."""
    rows = family.merge(contexts[['gene_id', 'level_a', 'level_b', 'subjects']], on='gene_id', validate='many_to_many')
    rows['test_id'] = pd.Series([requested_id(row.feature_id, row.level_a, row.level_b) for row in rows.itertuples(index=False)], index=rows.index, dtype=object)
    rows['n_requested_subjects'] = rows.subjects.map(lambda value: len(json.loads(value)))
    if rows.test_id.duplicated().any():
        raise ValueError('unique original event/contrast requests required')
    return rows


def validate_shard(table, expected, recipe, receipt, cohort, variant, index):
    if receipt['complete'] is not True or receipt['cohort'] != cohort or receipt['marker_variant'] != variant or receipt['shard_index'] != index or receipt['shard_count'] != recipe['shard_count'] or receipt['settings'] != recipe['settings'] or receipt['code_hashes'] != recipe['code_hashes'] or receipt['requested_tests'] != len(expected) or receipt['completed_tests'] != len(expected) or len(table) != len(expected):
        raise ValueError('all completed uniform scientific-recipe/count-family shards required')
    if table.test_id.duplicated().any() or set(table.test_id) != set(expected.test_id):
        raise ValueError('every original requested event/contrast identity required')
    boolean = table.converged.astype(str).str.lower()
    if not boolean.isin(['true', 'false']).all():
        raise ValueError('explicit numerical availability required')
    table = table.copy()
    table['converged'] = boolean.eq('true')
    cache_hits = table.fit_cache_hit.astype(str).str.lower()
    if not cache_hits.isin(['true', 'false']).all() or int(cache_hits.eq('true').sum()) != receipt['exact_fit_reuse'] or receipt['exact_fit_reuse'] + receipt['distinct_cache_evaluations'] != len(table):
        raise ValueError('exact reuse accounting must preserve every requested record')
    values = table[['p_value', 'raw_p_value', 'F_reference_p_value']].to_numpy(float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any() or not table.loc[~table.converged, ['p_value', 'raw_p_value', 'F_reference_p_value']].eq(1.).all().all() or not table.loc[~table.converged, 'statistic'].eq(0.).all():
        raise ValueError('failed requests must retain p1 and statistic zero')
    if not table.model_version.eq(recipe['settings']['model_version']).all() or not table.cohort.eq(cohort).all() or not table.marker_variant.eq(variant).all():
        raise ValueError('uniform model/cohort/marker recipe required')
    keys = ['test_id', 'feature_id', 'gene_id', 'event_type', 'level_a', 'level_b', 'n_requested_subjects']
    check = table[keys].merge(expected[keys], on='test_id', validate='one_to_one', suffixes=('_fit', '_expected'))
    if any(not check[f'{key}_fit'].eq(check[f'{key}_expected']).all() for key in keys[1:]):
        raise ValueError('original event labels and requested subjects changed')
    if not np.isfinite(table.statistic.to_numpy(float)).all() or table.statistic.lt(0).any() or not table.raw_p_value.eq(table.p_value).all():
        raise ValueError('finite nonnegative statistics and unchanged native p-values required')
    if table.loc[table.converged, 'n_subjects'].lt(4).any() or table.n_subjects.gt(table.n_requested_subjects).any() or not np.isfinite(table.loc[table.converged, 'effect_size'].to_numpy(float)).all() or receipt['usable'] != int(table.converged.sum()):
        raise ValueError('declared fit availability or effect differs from output')
    return table


def load_cohort(root, recipe, cohort, variant):
    family_path, context_path = [root / name for name in ('declared_family.tsv.gz', 'contexts.tsv.gz')]
    if file_hash(family_path) != recipe['family_sha256'] or file_hash(context_path) != recipe['contexts_sha256']:
        raise ValueError('original declared catalogue or subject contexts changed')
    family = pd.read_csv(family_path, sep='\t')
    contexts = pd.read_csv(context_path, sep='\t')
    contexts = contexts.loc[contexts.cohort.eq(cohort)]
    tables, receipts = [], []
    for index in range(recipe['shard_count']):
        folder = root / f'{variant}/{cohort}/shard_{index}'
        receipt_path, tests_path = folder / 'manifest.json', folder / 'tests.tsv.gz'
        receipt = json.loads(receipt_path.read_text())
        if receipt['recipe_sha256'] != file_hash(root / 'recipe.json') or receipt['tests_sha256'] != file_hash(tests_path):
            raise ValueError('original inference recipe and complete fit outputs required')
        expected = expected_requests(family.loc[family.shard_index.eq(index)], contexts.loc[contexts.shard_index.eq(index)])
        table = pd.read_csv(tests_path, sep='\t')
        tables.append(validate_shard(table, expected, recipe, receipt, cohort, variant, index))
        receipts.append(receipt)
    table = pd.concat(tables, ignore_index=True)
    if table.test_id.duplicated().any() or len(table) != recipe['requested_tests_by_cohort'][cohort]:
        raise ValueError('entire original cohort family required')
    return table, receipts


def split_table(table, method, tail, effect_column='effect_size'):
    local = normalize_pairs(table)
    local['method'] = method
    local['p_value'] = local[tail]
    local['raw_p_value'] = local[tail]
    local['coverage'] = local.median_gene_umis
    local['effect_size'] = local[effect_column]
    local['effect_vector'] = local.effect_size.map(lambda value: [float(value)])
    local['effect_features'] = [['inclusion'] for _ in range(len(local))]
    return local


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--marker-variant', choices=('all', 'junction'), default='all')
    parser.add_argument('--event-catalog', type=Path, required=True)
    parser.add_argument('--matrix-dir', type=Path, required=True)
    parser.add_argument('--gtf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError('new complete endpoint output required')
    recipe_path = args.root / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    if recipe['settings'] != FIT_SETTINGS or recipe['shard_count'] != 64 or recipe['declared_events'] != 23141:
        raise ValueError('complete original full-catalogue numerical recipe required')
    # Validate every cohort before creating an endpoint directory.
    cohorts, receipts = {}, {}
    for cohort in COHORTS:
        cohorts[cohort], receipts[cohort] = load_cohort(args.root, recipe, cohort, args.marker_variant)
    full = cohorts['full']
    catalogue = pd.read_csv(args.event_catalog, sep='\t')
    declared = pd.read_csv(args.root / 'declared_family.tsv.gz', sep='\t')
    if catalogue.feature_id.duplicated().any() or set(catalogue.feature_id) != set(declared.feature_id):
        raise ValueError('all original annotation event definitions required for LR mapping')
    input_paths = [args.event_catalog, args.gtf, *(args.matrix_dir / name for name in ('matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz'))]
    input_hashes = {str(path.resolve()): file_hash(path) for path in input_paths}
    method = f'Local read mixed binomial, {args.marker_variant} markers, integrated subject uncertainty'
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True)
    summaries = []
    for cohort, table in cohorts.items():
        table.to_csv(args.output_dir / f'{cohort}_tests.tsv.gz', sep='\t', index=False)
        censoring = table.effect_estimate_censored.astype(str).str.lower().eq('true')
        cache_hits = table.fit_cache_hit.astype(str).str.lower().eq('true')
        summaries.append(dict(cohort=cohort, requested=len(table), usable=int(table.converged.sum()), nominal_05=int(table.p_value.le(.05).sum()), summed_fit_runtime_seconds=float(table.runtime_seconds.sum()), exact_fit_reuse=int(cache_hits.sum()), distinct_cache_evaluations=len(table) - int(cache_hits.sum()), censoring=int(censoring.sum())))
    pd.DataFrame(summaries).to_csv(args.output_dir / 'fit_summary.tsv', sep='\t', index=False)
    correlations = []
    for label, tail in (('native_chi1', 'p_value'), ('F_reference', 'F_reference_p_value')):
        output = args.output_dir / label
        output.mkdir()
        name = method + (', F-reference sensitivity' if label == 'F_reference' else ', native chi1 reference')
        folds = [split_table(cohorts[cohort], name, tail) for cohort in ('fold0', 'fold1')]
        split_assessment(folds, repo, output, name)
        reporting_output = output / 'pooled_marker_reporting'
        reporting_output.mkdir()
        reporting_name = name + ', pooled marker reporting, unchanged tests'
        reporting_folds = [split_table(cohorts[cohort], reporting_name, tail, 'pooled_marker_effect') for cohort in ('fold0', 'fold1')]
        split_assessment(reporting_folds, repo, reporting_output, reporting_name)
        for fold, table in enumerate(folds):
            correlations.append(dict(tail=label, fold=fold, **coverage_correlation(table.p_value, table.coverage, table[['n_subjects']])))
    pd.DataFrame(correlations).to_csv(args.output_dir / 'coverage_correlations.tsv', sep='\t', index=False)
    full = full.assign(method=method)
    full_path = args.output_dir / 'full_tests.tsv.gz'
    full.to_csv(full_path, sep='\t', index=False)
    mapping_path = args.output_dir / 'lr_mapping.tsv.gz'
    command = [sys.executable, str(repo / 'extra_scripts/assess_event_tilgner_replication.py'), '--tests', str(full_path), '--events', str(args.event_catalog), '--tilgner-matrix', str(args.matrix_dir), '--gtf', str(args.gtf), '--output', str(mapping_path), '--summary', str(args.output_dir / 'lr_depth_summary.tsv'), '--top-per-contrast', '0']
    subprocess.run(command, check=True)
    mapped = pd.read_csv(mapping_path, sep='\t')
    # The same LR-evaluable identities/depth/direction policy for both tails.
    valid = mapped.mapping_complete.astype(str).str.lower().eq('true') & mapped.minimum_pooled_depth.ge(20) & mapped.pooled_replicated.notna()
    eligible = mapped.loc[valid].copy()
    eligible['pooled_replicated'] = eligible.pooled_replicated.astype(str).str.lower().eq('true')
    f_values = full.set_index(['feature_id', 'contrast_id']).F_reference_p_value
    for label in ('native_chi1', 'F_reference'):
        local = eligible.copy()
        if label == 'F_reference':
            local['p_value'] = [f_values.loc[(row.feature_id, row.contrast_id)] for row in local.itertuples(index=False)]
            local['raw_p_value'] = local.p_value
        local['method'] = method + (', F-reference sensitivity' if label == 'F_reference' else ', native chi1 reference')
        ranked = rank_direction_table(local, len(local))
        output = args.output_dir / label
        pd.DataFrame(ranked_direction_summary(ranked)).to_csv(output / 'lr_rank_summary.tsv', sep='\t', index=False)
        pd.DataFrame(ranked_category_summary(ranked)).to_csv(output / 'lr_event_type_composition.tsv', sep='\t', index=False)
        ranked.head(200).to_csv(output / 'lr_rank.tsv.gz', sep='\t', index=False)
        reported = reexpress_event_directions(local, full, 'pooled_marker_effect', local.method.iloc[0] + ', pooled marker reporting, unchanged tests' if len(local) else method + ', pooled marker reporting, unchanged tests')
        reporting_rank = rank_direction_table(reported, len(reported))
        if ranked[['feature_id', 'contrast_id', 'rank']].to_records(index=False).tolist() != reporting_rank[['feature_id', 'contrast_id', 'rank']].to_records(index=False).tolist():
            raise ValueError('reporting cannot change the LR-evaluable family or test ranks')
        reporting_output = output / 'pooled_marker_reporting'
        pd.DataFrame(ranked_direction_summary(reporting_rank)).to_csv(reporting_output / 'lr_rank_summary.tsv', sep='\t', index=False)
        pd.DataFrame(ranked_category_summary(reporting_rank)).to_csv(reporting_output / 'lr_event_type_composition.tsv', sep='\t', index=False)
        reporting_rank.head(200).to_csv(reporting_output / 'lr_rank.tsv.gz', sep='\t', index=False)
    code_paths = [Path(__file__), repo / 'extra_scripts/assess_event_tilgner_replication.py', Path(split_assessment.__code__.co_filename), Path(normalize_pairs.__code__.co_filename), Path(rank_direction_table.__code__.co_filename)]
    manifest = dict(recipe_sha256=file_hash(recipe_path), recipe=recipe, cohorts=receipts, LR_input_hashes=input_hashes, assessment_code_hashes={str(path.resolve()): file_hash(path) for path in code_paths}, LR='all complete tested associations freshly mapped, unchanged depth20 and nonzero finite direction requirement; no native-prefix restriction or significance filter', split='fixed published matched gene/pair universes, Simes event/pair aggregation and conjunction gene BH; all original unavailable fits retained atp1', estimand='common local-marker log-odds contrast; direction is not absolute RNA PSI or a pooled RNA-usage effect', reporting='separate equal-primer pooled-marker reporting, unchanged tests/gene-pair universes/LR association identities and ranks, missing or zero alternative effects count as LR nonagreement; split direction agreement assessed separately', calibration_limitation='native chi1 and F sensitivity are working-model diagnostics, not certified extreme-tail/joint-gene-FDR inference; numerical and toy checks do not establish biological validity', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
