#!/usr/bin/env python3
"""Whole annotation/LR diagnostic of event marker sources outside its classes."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_event_local_read_support import decode_event, file_hash
from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, stable_identifier
from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.local_marker_identifiability import transcript_marker_classes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--events', type=Path, required=True)
    parser.add_argument('--gtf', type=Path, required=True)
    parser.add_argument('--matrix-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError('new whole-catalogue marker-scope output required')
    recipe_path = args.source / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    if recipe.get('full_catalog') is not True or recipe.get('library_union') is not True or len(recipe['declared_events']) != 23141:
        raise ValueError('whole original marker catalogue required')
    catalogue = pd.read_csv(args.events, sep='\t').set_index('feature_id', verify_integrity=True)
    if set(catalogue.index) != set(recipe['declared_events']) or file_hash(args.events) != recipe['source_hashes'][str(args.events)] or file_hash(args.gtf) != recipe['source_hashes'][str(args.gtf)]:
        raise ValueError('unchanged source event definitions and annotation required')
    geometries = {row['event_id']: decode_event(row) for row in recipe['events']}
    annotation = read_gtf_exons(args.gtf)
    matrix, features, columns = read_tilgner_matrix(args.matrix_dir, args.gtf)
    selected_columns = columns.loc[columns.tealeaf_cell_type.notna() & columns.replicate.isin([1, 2]), 'column'].to_numpy(int)
    features['pooled_UMIs'] = np.asarray(matrix[:, selected_columns].sum(axis=1)).ravel()
    features['stable_transcript_id'] = features.transcript_id.map(lambda value: stable_identifier(value) if pd.notna(value) else None)
    abundance = features.dropna(subset=['stable_transcript_id']).groupby(['stable_gene_id', 'stable_transcript_id']).pooled_UMIs.sum().to_dict()
    records = []
    for key, event in catalogue.sort_index().iterrows():
        geometry = geometries.get(key)
        if geometry is None:
            for variant in ('all', 'junction'):
                records.append(dict(feature_id=key, gene_id=event.gene_id, event_type=event.event_type, marker_variant=variant, status='unavailable geometry'))
            continue
        chains = annotation[event.gene_id]['transcripts']
        included, excluded = set(event.included.split(',')), set(event.excluded.split(','))
        if not included or not excluded or included & excluded or not (included | excluded) <= set(chains):
            raise ValueError('original disjoint event transcript classes required')
        outside = set(chains) - included - excluded
        weights = {transcript: abundance.get((stable_identifier(event.gene_id), stable_identifier(transcript)), 0.) for transcript in chains}
        gene_total, event_total = sum(weights.values()), sum(weights[transcript] for transcript in included | excluded)
        for variant in ('all', 'junction'):
            classes = transcript_marker_classes(geometry, chains, variant == 'junction')
            sources = {transcript for transcript, flags in classes.items() if any(flags)}
            outside_sources = sources & outside
            outside_UMIs = sum(weights[transcript] for transcript in outside_sources)
            source_UMIs = sum(weights[transcript] for transcript in sources)
            record = dict(feature_id=key, gene_id=event.gene_id, event_type=event.event_type, marker_variant=variant, status='ok', annotated_transcripts=len(chains), included_transcripts=len(included), excluded_transcripts=len(excluded), outside_transcripts=len(outside), outside_marker_sources=len(outside_sources), outside_included_marker_sources=sum(classes[transcript][0] for transcript in outside), outside_excluded_marker_sources=sum(classes[transcript][1] for transcript in outside), outside_both_marker_sources=sum(all(classes[transcript]) for transcript in outside), total_LR_gene_UMIs=float(gene_total), total_LR_event_class_UMIs=float(event_total), outside_marker_source_LR_UMIs=float(outside_UMIs), outside_source_fraction_of_gene=float(outside_UMIs / gene_total) if gene_total else np.nan, outside_source_fraction_of_marker_source_RNA=float(outside_UMIs / source_UMIs) if source_UMIs else np.nan, event_class_fraction_of_gene=float(event_total / gene_total) if gene_total else np.nan)
            records.append(record)
    table = pd.DataFrame(records)
    if len(table) != 2 * len(recipe['declared_events']) or table.duplicated(['feature_id', 'marker_variant']).any():
        raise ValueError('every original event and marker variant required')
    summaries = []
    for scope, grouping in (('whole catalogue', ['marker_variant']), ('event type', ['marker_variant', 'event_type'])):
        for keys, frame in table.groupby(grouping):
            keys = keys if isinstance(keys, tuple) else (keys,)
            valid = frame.status.eq('ok')
            measurable = valid & frame.total_LR_gene_UMIs.ge(20)
            summaries.append(dict(zip(grouping, keys), scope=scope, declared_events=len(frame), geometries=int(valid.sum()), events_with_outside_marker_source=int((valid & frame.outside_marker_sources.gt(0)).sum()), events_with_gene_LR_depth_20=int(measurable.sum()), events_with_expressed_outside_marker_source=int((measurable & frame.outside_marker_source_LR_UMIs.gt(0)).sum()), median_outside_source_fraction_of_gene=float(frame.loc[measurable, 'outside_source_fraction_of_gene'].median()), median_outside_source_fraction_of_marker_source_RNA=float(frame.loc[measurable, 'outside_source_fraction_of_marker_source_RNA'].median()), median_event_class_fraction_of_gene=float(frame.loc[measurable, 'event_class_fraction_of_gene'].median())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / 'event_marker_scope.tsv.gz', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / 'summary.tsv', sep='\t', index=False)
    paths = [recipe_path, args.events, args.gtf, *(args.matrix_dir / name for name in ('matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz'))]
    manifest = dict(input_hashes={str(path): file_hash(path) for path in paths}, code_hashes={str(path): file_hash(path) for path in (Path(__file__), Path(transcript_marker_classes.__code__.co_filename))}, declared_events=len(recipe['declared_events']), scope='whole catalogue geometry/expression diagnostic, no p-value, LR-sign, discovery or rank selection', interpretation='outside-class annotated RNAs with potential marker sources, not estimated short-read contamination or protocol-specific capture probabilities; whole-chain marker flags do not mean one read spans every feature', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
