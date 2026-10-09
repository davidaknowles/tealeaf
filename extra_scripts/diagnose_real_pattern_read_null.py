#!/usr/bin/env python3
"""Diagnose a complete real-pattern null family without success filtering."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.audit_real_pattern_read_null import LAWS, trial_rng, validate_trials
from tealeaf.sc.local_read_diagnostics import local_read_design_diagnostic, local_read_failure_reason
from tealeaf.sc.local_read_null import local_read_pattern_null


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--null-root', type=Path, required=True)
    parser.add_argument('--assessment', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError('new diagnostic output required')
    paths = [args.null_root / 'recipe.json', args.null_root / 'parents.tsv.gz', args.assessment / 'manifest.json', args.assessment / 'tests.tsv.gz']
    recipe, receipt = [json.loads(path.read_text()) for path in (paths[0], paths[2])]
    if receipt['recipe'] != recipe or file_hash(paths[1]) != recipe['parents_sha256'] or len(receipt['receipts']) != recipe['shard_count'] or any(row['complete'] is not True for row in receipt['receipts']):
        raise ValueError('entire completed original null family required')
    parents = pd.read_csv(paths[1], sep='\t')
    frames = []
    for index, shard in enumerate(receipt['receipts']):
        path = args.null_root / f'shard_{index}/tests.tsv.gz'
        if shard['shard_index'] != index or shard['recipe_sha256'] != file_hash(paths[0]) or file_hash(path) != shard['tests_sha256']:
            raise ValueError('all original terminal shard identities and output hashes required')
        paths.append(path)
        frames.append(pd.read_csv(path, sep='\t'))
    table = validate_trials(pd.concat(frames, ignore_index=True), parents)
    if len(table) != recipe['requested_trials']:
        raise ValueError('all declared trials required, never successful fits only')
    parent_design, drawn_design = {}, {}
    for parent in parents.itertuples(index=False):
        original = np.asarray(json.loads(parent.counts_json), dtype=np.int64)
        parent_design[parent.parent_id] = local_read_design_diagnostic(original)
        for law_index, (law, sd) in enumerate(LAWS.items()):
            for draw in range(recipe['draws']):
                counts = local_read_pattern_null(original, trial_rng(parent.parent_id, law_index, draw), sd)
                drawn_design[parent.parent_id, law, draw] = local_read_design_diagnostic(counts)
    rows = []
    for row in table.itertuples(index=False):
        rows.append(dict(parent_id=row.parent_id, cohort=row.cohort, law=row.law, draw=row.draw, event_type=row.event_type, original_design_status=parent_design[row.parent_id]['design_status'], original_design_ready=parent_design[row.parent_id]['design_ready'], **drawn_design[row.parent_id, row.law, row.draw], failure_reason=local_read_failure_reason(row), converged=bool(row.converged), p_value=row.p_value, F_reference_p_value=row.F_reference_p_value))
    diagnosed = pd.DataFrame(rows)
    summaries = []
    for scope, keys in (('original design', ['cohort', 'law', 'original_design_status']), ('drawn design', ['cohort', 'law', 'design_status']), ('availability cause', ['cohort', 'law', 'failure_reason'])):
        for labels, frame in diagnosed.groupby(keys, sort=True):
            summaries.append(dict(zip(keys, labels), scope=scope, requested=len(frame), parents=frame.parent_id.nunique(), usable=int(frame.converged.sum()), native_05=int(frame.p_value.le(.05).sum()), native_01=int(frame.p_value.le(.01).sum()), F_05=int(frame.F_reference_p_value.le(.05).sum()), F_01=int(frame.F_reference_p_value.le(.01).sum())))
    args.output_dir.mkdir(parents=True)
    diagnosed.to_csv(args.output_dir / 'diagnosed_trials.tsv.gz', sep='\t', index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.output_dir / 'summary.tsv', sep='\t', index=False, na_rep='NA')
    code_paths = (Path(__file__), Path(local_read_design_diagnostic.__code__.co_filename), Path(validate_trials.__code__.co_filename), Path(local_read_pattern_null.__code__.co_filename))
    manifest = dict(input_hashes={str(path): file_hash(path) for path in paths}, code_hashes={str(path): file_hash(path) for path in code_paths}, trials=len(diagnosed), scope='complete unchanged null family, pre-draw original marker design and post-draw numerical diagnostics kept separate; no test/rank filtering or replacement endpoint', limitation='original marker eligibility is a descriptive fixed-parent stratum, numerical-success strata are selected outcomes and cannot certify calibration; repeated draws do not certify extreme tails or gene FDR', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary.loc[summary.scope.eq('original design')].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
