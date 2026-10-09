#!/usr/bin/env python3
"""Describe outside-class sources in unchanged published LR rank prefixes."""

import argparse
import json
from pathlib import Path

import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from tealeaf.sc.replication_scope import diagnose_ranked_marker_scope as diagnose


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scope', type=Path, required=True)
    parser.add_argument('--ranks', nargs=2, action='append', metavar=('LABEL', 'PATH'), required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError('new diagnostic output required')
    scope = pd.read_csv(args.scope, sep='\t')
    receipts, summaries, joined = {}, [], []
    for label, name in args.ranks:
        path = Path(name)
        ranked = pd.read_csv(path, sep='\t')
        summary, table = diagnose(ranked, scope, label)
        summaries.append(summary)
        joined.append(table)
        receipts[label] = dict(path=str(path), sha256=file_hash(path), original_ranked_associations=len(ranked))
    args.output_dir.mkdir(parents=True)
    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(args.output_dir / 'summary.tsv', sep='\t', index=False)
    pd.concat(joined, ignore_index=True).to_csv(args.output_dir / 'unchanged_prefix_scope.tsv.gz', sep='\t', index=False)
    manifest = dict(scope_sha256=file_hash(args.scope), inputs=receipts, code_hashes={str(path): file_hash(path) for path in (Path(__file__), Path(diagnose.__code__.co_filename))}, interpretation='descriptive annotation/RNA source diagnostic on unchanged complete/short method-own rank prefixes; no reranking, event exclusion, replacement LR endpoint, contamination estimation, causal attribution or production change', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary.loc[summary.cutoff.eq(100) & summary.selection.eq('all original prefix')].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
