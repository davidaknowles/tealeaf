#!/usr/bin/env python3
"""Whole marker catalogue inference, never a significance-selected pilot."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np
import pandas as pd
from scipy import stats

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs
from extra_scripts.collate_full_local_read_counts import KEYS
from tealeaf.sc.local_read_support import local_read_count_tensor
from tealeaf.sc import local_read_mixed


PRIMERS = ('poly(dT)', 'random hexamer')
COHORTS = ('fold0', 'fold1', 'full')
VARIANTS = ('all', 'junction')
FIT_SETTINGS = dict(model_version=local_read_mixed.MODEL_VERSION, min_gene_keys=25, min_subjects=4, max_iter=150, quadrature_tolerance=1e-3, node_schedule=[11, 21, 41], primary_tail='native chi-square 1, exploratory until real-count null and full-family checks', alternative_tail='F(1,measured subjects-1), diagnostic only', mean_prior='none')
FIELDS = ('test_id', 'feature_id', 'block_id', 'gene_id', 'event_type', 'effect', 'level_a', 'level_b', 'contrast_id', 'cohort', 'method', 'marker_variant', 'model_version', 'p_value', 'raw_p_value', 'F_reference_p_value', 'statistic', 'converged', 'n_subjects', 'n_requested_subjects', 'median_gene_umis', 'n_local_included_keys', 'n_local_excluded_keys', 'effect_size', 'effect_coordinate', 'counts_sha256', 'runtime_seconds', 'error', 'n_profiled_constant_primers', 'profiled_constant_primers', 'log_odds_effect', 'effect_estimate_censored', 'nuisance_parameter_boundary', 'null_baseline_sd', 'null_subject_sd', 'alternative_baseline_sd', 'alternative_subject_sd', 'null_objective', 'alternative_objective', 'quadrature_error', 'parameter_boundary', 'quadrature_orders_tried', 'final_quadrature_order', 'integration_policy')


def cohort_metadata(groups, folds, cohort):
    """Same fixed subject split and retained-cell groups in all variants."""
    if cohort not in COHORTS or groups.duplicated(['subject', 'cell_type', 'primer']).any() or folds.subject.duplicated().any() or set(folds.fold) != {0, 1}:
        raise ValueError('unique source groups and fixed two-fold subject declaration required')
    if not set(groups.subject) <= set(folds.subject) or set(groups.primer) != set(PRIMERS) or not groups.groupby(['subject', 'cell_type']).primer.apply(lambda values: set(values) == set(PRIMERS)).all():
        raise ValueError('every production pseudobulk needs its declared subject and both primers')
    metadata = groups[['subject', 'cell_type']].drop_duplicates().merge(folds[['subject', 'condition', 'fold']], on='subject', validate='many_to_one')
    if cohort != 'full':
        metadata = metadata.loc[metadata.fold.eq(int(cohort[-1]))].copy()
    metadata['mouse'] = metadata.subject
    return metadata.sort_values(['subject', 'cell_type']).reset_index(drop=True)


def gene_contexts(metadata, gene_counts):
    """Gene coverage only, no marker fit, event support or LR selection."""
    totals = gene_counts.groupby(['subject', 'cell_type']).gene_keys.sum()
    coverage = np.array([totals.get((row.subject, row.cell_type), 0.) for row in metadata.itertuples(index=False)])
    contexts = []
    for specification, _ in covered_celltype_pairwise_designs(metadata, coverage, min_gene_umis=FIT_SETTINGS['min_gene_keys'], min_samples=0, min_celltype_mice=FIT_SETTINGS['min_subjects']):
        rows, local, _, _, _, levels = specification
        contexts.append(dict(level_a=levels[0], level_b=levels[1], subjects=json.dumps(sorted(local.mouse.astype(str).unique())), median_gene_umis=float(np.median(coverage[rows]))))
    return contexts


def scientific_code_hashes():
    # Capture code before work, not after a potentially long worker has run.
    paths = (Path(__file__), Path(local_read_mixed.__file__), Path(local_read_count_tensor.__code__.co_filename), Path(covered_celltype_pairwise_designs.__code__.co_filename))
    return {str(path.resolve()): file_hash(path) for path in paths}


def prepare(args):
    if args.output_root.exists():
        raise ValueError('new whole-family inference recipe required')
    source_paths = [args.source / name for name in ('manifest.json', 'groups.tsv', 'declared_family.tsv.gz', 'gene_keys.tsv.gz', 'feature_definitions.tsv')]
    manifest = json.loads(source_paths[0].read_text())
    shards = manifest['count_shards']
    if len(manifest['library_receipts']) != 8 or any(receipt['complete'] is not True for receipt in manifest['library_receipts']) or manifest['production_cell_qc']['exact_cached_group_and_primer_total_match'] is not True or [row['shard_index'] for row in shards] != list(range(len(shards))):
        raise ValueError('whole frozen library/count-shard catalogue required')
    groups = pd.read_csv(source_paths[1], sep='\t', dtype={'subject': str})
    family = pd.read_csv(source_paths[2], sep='\t')
    genes = pd.read_csv(source_paths[3], sep='\t', dtype={'subject': str})
    folds = pd.read_csv(args.subject_folds, sep='\t', dtype={'subject': str})
    if family.feature_id.duplicated().any() or len(family) != manifest['declared_events'] or not set(genes.gene_id) <= set(family.gene_id):
        raise ValueError('complete declared event family and gene molecule counts required')
    metadata = {cohort: cohort_metadata(groups, folds, cohort) for cohort in COHORTS}
    gene_tables = dict(tuple(genes.groupby('gene_id', sort=False)))
    empty = genes.iloc[:0]
    contexts = []
    for gene, local in family.groupby('gene_id', sort=True):
        if local.shard_index.nunique() != 1:
            raise ValueError('a gene must remain within one declared count shard')
        for cohort in COHORTS:
            for context in gene_contexts(metadata[cohort], gene_tables.get(gene, empty)):
                contexts.append(dict(gene_id=gene, cohort=cohort, shard_index=int(local.shard_index.iloc[0]), n_events=len(local), **context))
    contexts = pd.DataFrame(contexts, columns=['gene_id', 'cohort', 'shard_index', 'n_events', 'level_a', 'level_b', 'subjects', 'median_gene_umis'])
    args.output_root.mkdir(parents=True)
    for path in (source_paths[2], args.subject_folds):
        shutil.copyfile(path, args.output_root / path.name)
    contexts.to_csv(args.output_root / 'contexts.tsv.gz', sep='\t', index=False)
    requests = {cohort: int(contexts.loc[contexts.cohort.eq(cohort), 'n_events'].sum()) for cohort in COHORTS}
    receipt = dict(source=str(args.source.resolve()), input_hashes={str(path.resolve()): file_hash(path) for path in (*source_paths, args.subject_folds)}, code_hashes=scientific_code_hashes(), settings=FIT_SETTINGS, declared_events=len(family), shard_count=len(shards), count_shards=shards, requested_tests_by_cohort=requests, contexts_sha256=file_hash(args.output_root / 'contexts.tsv.gz'), family_sha256=file_hash(args.output_root / 'declared_family.tsv.gz'), scope='whole-catalogue inference preparation, gene coverage only; all event/contrast failures remain p1, no significance or LR outcome selection', production_changes=False)
    (args.output_root / 'recipe.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(declared_events=len(family), requested_tests_by_cohort=requests, contexts=len(contexts)), indent=2), flush=True)


def requested_id(feature, first, second):
    return f'{feature}|cell_type|{first}|{second}'


def fit_record(event, context, lookup, cohort, variant):
    subjects = json.loads(context.subjects)
    values = local_read_count_tensor(event.feature_id, subjects, (context.level_a, context.level_b), PRIMERS, lookup)
    record = dict(test_id=requested_id(event.feature_id, context.level_a, context.level_b), feature_id=event.feature_id, block_id=event.feature_id.removeprefix('SUPPA2:'), gene_id=event.gene_id, event_type=event.event_type, effect='cell_type', level_a=context.level_a, level_b=context.level_b, contrast_id=f'cell_type__{context.level_a}__{context.level_b}', cohort=cohort, method=f'Local read mixed binomial, {variant} markers', marker_variant=variant, model_version=local_read_mixed.MODEL_VERSION, p_value=1., raw_p_value=1., F_reference_p_value=1., statistic=0., converged=False, n_subjects=0, n_requested_subjects=len(subjects), median_gene_umis=context.median_gene_umis, n_local_included_keys=int(values[..., 0].sum()), n_local_excluded_keys=int(values[..., 1].sum()), effect_size=np.nan, effect_coordinate='local marker inclusion log-odds, not RNA PSI', counts_sha256=hashlib.sha256(values.tobytes()).hexdigest(), error='')
    start = time.monotonic()
    try:
        result = local_read_mixed.local_read_mixed_adaptive_test(local_read_mixed.LocalReadMixed(values), node_schedule=tuple(FIT_SETTINGS['node_schedule']), max_iter=FIT_SETTINGS['max_iter'], quadrature_tolerance=FIT_SETTINGS['quadrature_tolerance'])
        record.update({key: value for key, value in result.items() if key in FIELDS})
        record['raw_p_value'] = record['p_value']
        if result['converged']:
            record['effect_size'] = result['log_odds_effect']
            record['F_reference_p_value'] = float(stats.f.sf(result['statistic'], 1, result['n_subjects'] - 1))
    except (ValueError, np.linalg.LinAlgError) as exc:
        record['error'] = str(exc)
    record['runtime_seconds'] = time.monotonic() - start
    return record


def fit(args):
    recipe_path = args.output_root / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    if recipe['settings'] != FIT_SETTINGS or recipe['code_hashes'] != scientific_code_hashes() or args.cohort not in COHORTS or args.marker_variant not in VARIANTS or not 0 <= args.shard_index < recipe['shard_count']:
        raise ValueError('declared uniform scientific recipe/cohort/count shard required')
    family_path, contexts_path = [args.output_root / name for name in ('declared_family.tsv.gz', 'contexts.tsv.gz')]
    if file_hash(family_path) != recipe['family_sha256'] or file_hash(contexts_path) != recipe['contexts_sha256']:
        raise ValueError('original declared family/subject contexts changed')
    source = Path(recipe['source'])
    manifest_path = source / 'manifest.json'
    if file_hash(manifest_path) != recipe['input_hashes'][str(manifest_path.resolve())]:
        raise ValueError('full count receipt changed')
    shard = recipe['count_shards'][args.shard_index]
    counts_path = source / f'shard_{args.shard_index}/counts.tsv.gz'
    if file_hash(counts_path) != shard['counts_sha256']:
        raise ValueError('frozen local counts changed')
    folder = args.output_root / f'{args.marker_variant}/{args.cohort}/shard_{args.shard_index}'
    if folder.exists():
        raise ValueError('never overwrite a completed or partial whole-family worker')
    family = pd.read_csv(family_path, sep='\t')
    family = family.loc[family.shard_index.eq(args.shard_index)]
    contexts = pd.read_csv(contexts_path, sep='\t')
    contexts = contexts.loc[contexts.shard_index.eq(args.shard_index) & contexts.cohort.eq(args.cohort)]
    counts = pd.read_csv(counts_path, sep='\t', dtype={'subject': str})
    if len(counts) != shard['count_rows'] or not set(counts.feature_id) <= set(family.feature_id) or counts.duplicated(KEYS).any():
        raise ValueError('unique compatible local count-shard rows required')
    fields = ['included', 'excluded'] if args.marker_variant == 'all' else ['junction_included', 'junction_excluded']
    lookup = counts[KEYS + fields].rename(columns=dict(zip(fields, ['included', 'excluded']))).set_index(KEYS).to_dict('index')
    events = {gene: local for gene, local in family.groupby('gene_id', sort=False)}
    requested = int(contexts.n_events.sum())
    folder.mkdir(parents=True)
    completed, usable, start = 0, 0, time.monotonic()
    # Stream rows for inspection, but only a terminal receipt permits collation.
    with gzip.open(folder / 'tests.tsv.gz', 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        for context in contexts.itertuples(index=False):
            local = events[context.gene_id]
            if len(local) != context.n_events:
                raise ValueError('declared event family differs from its context')
            for event in local.itertuples(index=False):
                record = fit_record(event, context, lookup, args.cohort, args.marker_variant)
                writer.writerow(record)
                completed += 1
                usable += int(record['converged'])
                if completed % 100 == 0:
                    handle.flush()
                    print(f'{completed}/{requested} requests, {usable} usable, {time.monotonic() - start:.1f} seconds', flush=True)
    if completed != requested:
        raise ValueError('every declared event/contrast request must be retained')
    receipt = dict(recipe_sha256=file_hash(recipe_path), code_hashes=recipe['code_hashes'], settings=FIT_SETTINGS, cohort=args.cohort, marker_variant=args.marker_variant, shard_index=args.shard_index, shard_count=recipe['shard_count'], declared_events=len(family), requested_tests=requested, completed_tests=completed, usable=usable, elapsed_seconds=time.monotonic() - start, tests_sha256=file_hash(folder / 'tests.tsv.gz'), complete=True, production_changes=False)
    (folder / 'manifest.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'fit'))
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--subject-folds', type=Path)
    parser.add_argument('--cohort', choices=COHORTS)
    parser.add_argument('--marker-variant', choices=VARIANTS, default='all')
    parser.add_argument('--shard-index', type=int)
    args = parser.parse_args()
    if args.mode == 'prepare':
        if args.source is None or args.subject_folds is None:
            parser.error('complete count family and frozen subject split required')
        prepare(args)
    else:
        if args.cohort is None or args.shard_index is None:
            parser.error('declared cohort and shard required')
        fit(args)


if __name__ == '__main__':
    main()
