#!/usr/bin/env python3
"""Matched numerical-boundary audit, not a selected-panel power benchmark."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.assess_library_read_model_pilot import validate_family
from extra_scripts.audit_event_local_read_support import file_hash


KEYS = ['fold', 'test_id', 'model', 'variant']


def matched_pilots(old, new, cases, versions=('v1', 'v2')):
    validate_family(old, cases)
    validate_family(new, cases)
    table = old.merge(new, on=KEYS, suffixes=('_old', '_new'), validate='one_to_one')
    if not table.counts_sha256_old.eq(table.counts_sha256_new).all() or not table.requested_subjects_old.eq(table.requested_subjects_new).all():
        raise ValueError('original counts and requested subjects must remain identical')
    if len(versions) != 2 or any(version not in ('v1', 'v2') for version in versions):
        raise ValueError('two explicit supported model versions required')
    for frame, version in zip((old, new), versions):
        observed = frame.loc[frame.model.eq('unconditional'), 'model_version'].dropna().unique()
        if set(observed) != {f'local_read_binomial_random_intercept_slope_{version}'}:
            raise ValueError('explicit original and updated unconditional model versions required')
    return table


def failure_reason(row):
    if row.converged:
        return 'usable'
    if pd.notna(row.error) and str(row.error):
        return str(row.error)
    tolerance = 1e-4 if row.model == 'conditional' else 1e-3
    if getattr(row, 'quadrature_error', 0.) > tolerance:
        return 'doubled-order quadrature mismatch'
    if getattr(row, 'parameter_boundary', False) is True or str(getattr(row, 'parameter_boundary', '')).lower() == 'true':
        return 'nonvalidated parameter boundary'
    return 'other numerical unavailability'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old', type=Path, required=True)
    parser.add_argument('--new', type=Path, required=True)
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--old-version', choices=('v1', 'v2'), default='v1')
    parser.add_argument('--new-version', choices=('v1', 'v2'), default='v2')
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError('new paired audit output required')
    old, new = [pd.read_csv(folder / 'tests.tsv.gz', sep='\t') for folder in (args.old, args.new)]
    receipts = [json.loads((folder / 'manifest.json').read_text()) for folder in (args.old, args.new)]
    if receipts[0]['source_hashes'] != receipts[1]['source_hashes'] or len(receipts[0]['shards']) != 16 or len(receipts[1]['shards']) != 16:
        raise ValueError('two complete pilots on identical source counts required')
    cases = pd.read_csv(args.cases, sep='\t')
    table = matched_pilots(old, new, cases, (args.old_version, args.new_version))
    summaries = []
    groups = ['fold', 'panel_new', 'model', 'variant']
    for keys, local in table.groupby(groups):
        a, b = local.converged_old, local.converged_new
        both = a & b
        finite = both & np.isfinite(local.log_odds_effect_old) & np.isfinite(local.log_odds_effect_new)
        summaries.append(dict(zip(['fold', 'panel', 'model', 'variant'], keys), requested=len(local), old_usable=int(a.sum()), new_usable=int(b.sum()), recovered=int((~a & b).sum()), newly_unavailable=int((a & ~b).sum()), same_available_direction_agreement=int((np.sign(local.loc[finite, 'log_odds_effect_old']) == np.sign(local.loc[finite, 'log_odds_effect_new'])).sum()), same_available_direction_denominator=int(finite.sum()), old_nominal_05=int(local.p_value_old.le(.05).sum()), new_nominal_05=int(local.p_value_new.le(.05).sum())))
    failures = []
    for name, frame in (('old', old), ('new', new)):
        local = frame.copy()
        local['failure_reason'] = [failure_reason(row) for row in local.itertuples(index=False)]
        detail = local.groupby(['fold', 'panel', 'model', 'variant', 'failure_reason']).size().rename('requested').reset_index()
        failures.append(detail.assign(version=name))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / 'matched_tests.tsv.gz', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / 'summary.tsv', sep='\t', index=False)
    pd.concat(failures).to_csv(args.output_dir / 'failure_reasons.tsv', sep='\t', index=False)
    paths = [args.cases, *(folder / name for folder in (args.old, args.new) for name in ('tests.tsv.gz', 'manifest.json'))]
    manifest = dict(input_hashes={str(path): file_hash(path) for path in paths}, requested_fits=len(table), declared_model_versions=[args.old_version, args.new_version], scope='matched numerical availability and direction audit of the unchanged selected panel; nominal calls are descriptive, not discoveries, calibrated power, unbiased split replication or own-ranked LR agreement', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
