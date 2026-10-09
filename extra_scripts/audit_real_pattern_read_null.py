#!/usr/bin/env python3
"""Frozen real-depth stress nulls for whole-catalogue local-marker inference."""

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.full_local_read_models import COHORTS, FIT_SETTINGS, PRIMERS, fit_record, requested_id, scientific_code_hashes
from tealeaf.sc.local_read_support import local_read_count_tensor
from tealeaf.sc.local_read_null import local_read_pattern_null, PATTERN_NULL_VERSION


SEED = 20261009
LAWS = dict(strict=0., biological=.8)
SHARDS = 32
DRAWS = 4


def stable_hash(value):
    return hashlib.sha256(f'{SEED}|{value}'.encode()).hexdigest()


def trial_rng(parent_id, law_index, draw):
    return np.random.default_rng(np.random.SeedSequence([SEED, int(stable_hash(parent_id)[:16], 16), law_index, draw]))


def choose_parents(family, contexts, cohort='full', per_stratum=16):
    """One frozen context/event per gene, no counts, fits or LR outcomes."""
    if family.feature_id.duplicated().any():
        raise ValueError('unique frozen source catalogue required')
    if cohort not in COHORTS:
        raise ValueError('declared original cohort required')
    local = contexts.loc[contexts.cohort.eq(cohort)].copy()
    local['context_id'] = [requested_id(row.gene_id, row.level_a, row.level_b) for row in local.itertuples(index=False)]
    if local.context_id.duplicated().any():
        raise ValueError('unique frozen gene/contrast contexts required')
    local['_choice'] = local.context_id.map(stable_hash)
    local = local.sort_values('_choice').drop_duplicates('gene_id').drop(columns='_choice')
    events = family.assign(_choice=family.feature_id.map(stable_hash)).sort_values('_choice').drop_duplicates('gene_id').drop(columns='_choice')
    local = local.merge(events[['gene_id', 'feature_id', 'event_type']], on='gene_id', validate='one_to_one')
    local['n_requested_subjects'] = local.subjects.map(lambda value: len(json.loads(value)))
    local['subject_stratum'] = np.where(local.n_requested_subjects.lt(8), '4-7', '8+')
    # Coverage gates/strata are pre-fit gene totals, not informative-class depth.
    local['coverage_quartile'] = pd.qcut(local.median_gene_umis.rank(method='first'), 4, labels=False)
    local['_choice'] = local.gene_id.map(stable_hash)
    local = local.sort_values('_choice').groupby(['subject_stratum', 'coverage_quartile'], sort=True).head(per_stratum).drop(columns='_choice').sort_values(['subject_stratum', 'coverage_quartile', 'gene_id']).reset_index(drop=True)
    local['parent_id'] = [f'{cohort}|{requested_id(row.feature_id, row.level_a, row.level_b)}' for row in local.itertuples(index=False)]
    return local


def prepare(args):
    if args.output_root.exists():
        raise ValueError('new null recipe required')
    recipe_path = args.root / 'recipe.json'
    source_recipe = json.loads(recipe_path.read_text())
    if source_recipe['settings'] != FIT_SETTINGS or source_recipe['code_hashes'] != scientific_code_hashes() or source_recipe['declared_events'] != 23141 or source_recipe['shard_count'] != 64:
        raise ValueError('original whole-catalogue inference preparation required')
    family_path, context_path = [args.root / name for name in ('declared_family.tsv.gz', 'contexts.tsv.gz')]
    if file_hash(family_path) != source_recipe['family_sha256'] or file_hash(context_path) != source_recipe['contexts_sha256']:
        raise ValueError('frozen source family/context hashes required')
    family, contexts = [pd.read_csv(path, sep='\t') for path in (family_path, context_path)]
    parents = pd.concat([choose_parents(family, contexts, cohort) for cohort in COHORTS], ignore_index=True)
    if len(parents) < SHARDS:
        raise ValueError('enough declared parents to populate every null shard required')
    chunks = []
    source = Path(source_recipe['source'])
    for index, packet in enumerate(source_recipe['count_shards']):
        path = source / f'shard_{index}/counts.tsv.gz'
        if file_hash(path) != packet['counts_sha256']:
            raise ValueError('original full count shard changed')
        frame = pd.read_csv(path, sep='\t', dtype={'subject': str})
        chunks.append(frame.loc[frame.feature_id.isin(parents.feature_id), ['feature_id', 'subject', 'cell_type', 'primer', 'included', 'excluded']])
    frame = pd.concat(chunks, ignore_index=True)
    lookup = frame.set_index(['feature_id', 'subject', 'cell_type', 'primer']).to_dict('index')
    tensors = [local_read_count_tensor(row.feature_id, json.loads(row.subjects), (row.level_a, row.level_b), PRIMERS, lookup).astype(np.int64) for row in parents.itertuples(index=False)]
    parents['counts_json'] = [json.dumps(values.tolist(), separators=(',', ':')) for values in tensors]
    parents['original_counts_sha256'] = [hashlib.sha256(values.tobytes()).hexdigest() for values in tensors]
    args.output_root.mkdir(parents=True)
    parents.to_csv(args.output_root / 'parents.tsv.gz', sep='\t', index=False)
    source_paths = [recipe_path, family_path, context_path]
    receipt = dict(source_root=str(args.root.resolve()), source_hashes={str(path.resolve()): file_hash(path) for path in source_paths}, audit_code_sha256=file_hash(Path(__file__)), fitting_code_hashes=scientific_code_hashes(), generator_version=PATTERN_NULL_VERSION, generator_sha256=file_hash(Path(local_read_pattern_null.__code__.co_filename)), settings=FIT_SETTINGS, parent_count=len(parents), draws=DRAWS, laws=LAWS, shard_count=SHARDS, seed=SEED, parents_sha256=file_hash(args.output_root / 'parents.tsv.gz'), requested_trials=len(parents) * DRAWS * len(LAWS), selection='one event and context per gene/cohort by frozen hash, up to16 genes per original subject-count/gene-coverage-quartile stratum in each cohort; all selected markerless and failed trials retained', law='fixed subject/primer baseline pooled across cell types, original per-type depths, fresh binomials and independent mean-zero shared-primer normal slopes; arbitrary baselines stress random-intercept misspecification, not that model own-law null or RNA PSI', production_changes=False)
    (args.output_root / 'recipe.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(parents=len(parents), trials=receipt['requested_trials']), indent=2), flush=True)


def fit(args):
    recipe_path = args.output_root / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    parent_path = args.output_root / 'parents.tsv.gz'
    if recipe['settings'] != FIT_SETTINGS or recipe['audit_code_sha256'] != file_hash(Path(__file__)) or recipe['fitting_code_hashes'] != scientific_code_hashes() or recipe['generator_version'] != PATTERN_NULL_VERSION or recipe['generator_sha256'] != file_hash(Path(local_read_pattern_null.__code__.co_filename)) or recipe['parents_sha256'] != file_hash(parent_path):
        raise ValueError('frozen full scientific null recipe required')
    if not 0 <= args.shard_index < SHARDS or recipe['shard_count'] != SHARDS:
        raise ValueError('declared null shard required')
    folder = args.output_root / f'shard_{args.shard_index}'
    if folder.exists():
        raise ValueError('never overwrite a complete or partial null shard')
    parents = pd.read_csv(parent_path, sep='\t').iloc[args.shard_index::SHARDS]
    rows = []
    for parent in parents.itertuples(index=False):
        original = np.asarray(json.loads(parent.counts_json), dtype=np.int64)
        if hashlib.sha256(original.tobytes()).hexdigest() != parent.original_counts_sha256:
            raise ValueError('original parent counts changed')
        subjects = json.loads(parent.subjects)
        context = SimpleNamespace(subjects=parent.subjects, level_a=parent.level_a, level_b=parent.level_b, median_gene_umis=parent.median_gene_umis)
        for law_index, (law, sd) in enumerate(LAWS.items()):
            for draw in range(DRAWS):
                observed = local_read_pattern_null(original, trial_rng(parent.parent_id, law_index, draw), sd)
                lookup = {(parent.feature_id, subject, level, primer): dict(included=int(observed[u, p, c, 0]), excluded=int(observed[u, p, c, 1])) for u, subject in enumerate(subjects) for p, primer in enumerate(PRIMERS) for c, level in enumerate((parent.level_a, parent.level_b))}
                record = fit_record(parent, context, lookup, parent.cohort, 'all')
                record.update(parent_id=parent.parent_id, law=law, draw=draw, generator_version=PATTERN_NULL_VERSION, original_counts_sha256=parent.original_counts_sha256, subject_stratum=parent.subject_stratum, coverage_quartile=parent.coverage_quartile)
                rows.append(record)
        print(f'{len(rows)} requested trials completed', flush=True)
    folder.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(folder / 'tests.tsv.gz', sep='\t', index=False)
    receipt = dict(recipe_sha256=file_hash(recipe_path), parents_sha256=recipe['parents_sha256'], shard_index=args.shard_index, parent_ids=parents.parent_id.tolist(), requested_trials=len(parents) * len(LAWS) * DRAWS, completed_trials=len(rows), tests_sha256=file_hash(folder / 'tests.tsv.gz'), complete=True, production_changes=False)
    (folder / 'manifest.json').write_text(json.dumps(receipt, indent=2) + '\n')


def validate_trials(frame, parents):
    expected = {(parent, law, draw) for parent in parents.parent_id for law in LAWS for draw in range(DRAWS)}
    if frame.duplicated(['parent_id', 'law', 'draw']).any() or set(frame[['parent_id', 'law', 'draw']].itertuples(index=False, name=None)) != expected:
        raise ValueError('every original null trial identity required')
    frame = frame.copy()
    frame['converged'] = frame.converged.astype(str).str.lower().map({'true': True, 'false': False})
    values = frame[['p_value', 'F_reference_p_value']].to_numpy(float)
    if frame.converged.isna().any() or not np.isfinite(values).all() or ((values < 0) | (values > 1)).any() or not frame.loc[~frame.converged, ['p_value', 'F_reference_p_value']].eq(1.).all().all() or not frame.loc[~frame.converged, 'statistic'].eq(0.).all():
        raise ValueError('failed trials must remain in denominators atp1')
    columns = ['feature_id', 'gene_id', 'cohort', 'level_a', 'level_b', 'original_counts_sha256', 'subject_stratum', 'coverage_quartile', 'n_requested_subjects']
    check = frame[['parent_id', *columns]].merge(parents[['parent_id', *columns]], on='parent_id', validate='many_to_one', suffixes=('_fit', '_original'))
    if any(not check[f'{column}_fit'].eq(check[f'{column}_original']).all() for column in columns) or not frame.model_version.eq(FIT_SETTINGS['model_version']).all() or not frame.generator_version.eq(PATTERN_NULL_VERSION).all():
        raise ValueError('frozen parent labels/counts/cohort and generating recipe required')
    expected_hashes = {}
    for parent in parents.itertuples(index=False):
        original = np.asarray(json.loads(parent.counts_json), dtype=np.int64)
        if hashlib.sha256(original.tobytes()).hexdigest() != parent.original_counts_sha256 or original.shape != (parent.n_requested_subjects, len(PRIMERS), 2, 2):
            raise ValueError('exact original subject/primer count tensor required')
        for law_index, (law, sd) in enumerate(LAWS.items()):
            for draw in range(DRAWS):
                observed = local_read_pattern_null(original, trial_rng(parent.parent_id, law_index, draw), sd)
                expected_hashes[(parent.parent_id, law, draw)] = hashlib.sha256(observed.tobytes()).hexdigest()
    if any(row.counts_sha256 != expected_hashes[(row.parent_id, row.law, row.draw)] for row in frame.itertuples(index=False)):
        raise ValueError('every actual count draw must match its frozen seed/parent/generator')
    return frame


def assess(args):
    if args.output_dir.exists():
        raise ValueError('new complete null assessment required')
    recipe_path = args.output_root / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    parent_path = args.output_root / 'parents.tsv.gz'
    if file_hash(parent_path) != recipe['parents_sha256'] or recipe['settings'] != FIT_SETTINGS or recipe['audit_code_sha256'] != file_hash(Path(__file__)) or recipe['generator_version'] != PATTERN_NULL_VERSION or recipe['draws'] != DRAWS or recipe['laws'] != LAWS or recipe['shard_count'] != SHARDS:
        raise ValueError('all original declared parents required')
    parents = pd.read_csv(parent_path, sep='\t')
    frames, receipts = [], []
    for index in range(SHARDS):
        folder = args.output_root / f'shard_{index}'
        receipt = json.loads((folder / 'manifest.json').read_text())
        local_parents = parents.iloc[index::SHARDS]
        expected_parents = local_parents.parent_id.tolist()
        path = folder / 'tests.tsv.gz'
        if receipt['complete'] is not True or receipt['recipe_sha256'] != file_hash(recipe_path) or receipt['parents_sha256'] != recipe['parents_sha256'] or receipt['shard_index'] != index or receipt['parent_ids'] != expected_parents or receipt['tests_sha256'] != file_hash(path) or receipt['requested_trials'] != len(expected_parents) * DRAWS * len(LAWS) or receipt['completed_trials'] != receipt['requested_trials']:
            raise ValueError('every complete frozen null shard required')
        frame = pd.read_csv(path, sep='\t')
        frame = validate_trials(frame, local_parents)
        frames.append(frame)
        receipts.append(receipt)
    table = pd.concat(frames, ignore_index=True)
    if len(table) != recipe['requested_trials']:
        raise ValueError('whole prespecified parent/trial family required')
    summaries = []
    for scope, grouping in (('all requested', ['cohort', 'law']), ('subject/coverage strata', ['cohort', 'law', 'subject_stratum', 'coverage_quartile'])):
        for keys, frame in table.groupby(grouping):
            keys = keys if isinstance(keys, tuple) else (keys,)
            for tail in ('p_value', 'F_reference_p_value'):
                for threshold in (.05, .01):
                    rejected = int(frame[tail].le(threshold).sum())
                    interval = binomtest(rejected, len(frame)).proportion_ci()
                    summaries.append(dict(zip(grouping, keys), scope=scope, tail=tail, threshold=threshold, requested=len(frame), usable=int(frame.converged.sum()), n_parents=frame.parent_id.nunique(), rejected=rejected, rejection_rate=rejected / len(frame), ci_low=interval.low, ci_high=interval.high))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / 'tests.tsv.gz', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / 'summary.tsv', sep='\t', index=False)
    manifest = dict(recipe=recipe, receipts=receipts, scope='all frozen real-depth pattern stress trials; repeated draws and shared design parents do not certify extreme tails or whole-gene FDR; failures retained atp1', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summaries).loc[lambda frame: frame.scope.eq('all requested')].to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'fit', 'assess'))
    parser.add_argument('--root', type=Path)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--shard-index', type=int)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    if args.mode == 'prepare' and args.root is None or args.mode == 'fit' and args.shard_index is None or args.mode == 'assess' and args.output_dir is None:
        parser.error('source root, declared shard or assessment output required for selected mode')
    globals()[args.mode](args)


if __name__ == '__main__':
    main()
