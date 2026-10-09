#!/usr/bin/env python3
"""Guard and partition a whole annotation-catalogue local-read count family."""

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash


KEYS = ["feature_id", "subject", "cell_type", "primer"]
COUNT_COLUMNS = ["included", "excluded", "junction_included", "junction_excluded", "gene_keys"]


def summarize_markers(signatures):
    """Unanimous class keys, with conflicts excluded and gene keys unchanged."""
    if signatures.duplicated([*KEYS, "signature"]).any():
        raise ValueError("unique aggregated signature identities required")
    values = signatures[["signature", "count"]].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or not (values == np.floor(values)).all() or (signatures.signature > 31).any():
        raise ValueError("exact nonnegative five-bit signatures and key counts required")
    frame = signatures[KEYS].copy()
    bits, counts = signatures.signature.to_numpy(dtype=int), signatures['count'].to_numpy(dtype=np.int64)
    for first, second, inc, exc in ((5, 10, "included", "excluded"), (4, 8, "junction_included", "junction_excluded")):
        a, b = (bits & first) != 0, (bits & second) != 0
        frame[inc], frame[exc] = counts * (a & ~b), counts * (b & ~a)
    frame['gene_keys'] = counts
    return frame.groupby(KEYS, sort=False)[COUNT_COLUMNS].sum().reset_index()


def gene_key_counts(markers, family):
    """Verify gene-exonic molecule totals agree across every event definition.

    Marker signatures depend on the local event, but the gene eligibility gate
    and exact CB/UB/gene key do not. Zero-count genes are kept in the family.
    """
    if family.feature_id.duplicated().any() or not set(markers.feature_id) <= set(family.feature_id):
        raise ValueError("declared unique annotation family required")
    lookup = family.set_index("feature_id").gene_id
    values = markers.assign(gene_id=markers.feature_id.map(lookup))
    group = ["gene_id", "subject", "cell_type", "primer"]
    checks = values.groupby(group).gene_keys.agg(['min', 'max', 'size'])
    expected = family.groupby('gene_id').size()
    if not checks['min'].eq(checks['max']).all() or not checks['size'].eq(checks.index.get_level_values('gene_id').map(expected)).all():
        raise ValueError("gene-key accounting must be identical across every available event geometry")
    return checks['max'].rename('gene_keys').reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--shard-count', type=int, default=64)
    args = parser.parse_args()
    if args.output_dir.exists() or args.shard_count < 1:
        raise ValueError('new output and positive shard count required')
    recipe_path = args.source / 'recipe.json'
    recipe = json.loads(recipe_path.read_text())
    if not recipe.get('full_catalog') or not recipe.get('library_union') or len(recipe['bams']) != 8 or recipe['production_cell_qc']['exact_cached_group_and_primer_total_match'] is not True:
        raise ValueError('complete full catalogue and exact production cell scope required')
    recipe_hash = file_hash(recipe_path)
    declared_path, definitions_path = [args.source / name for name in ('declared_family.tsv.gz', 'feature_definitions.tsv')]
    family = pd.read_csv(declared_path, sep='\t')
    if family.feature_id.duplicated().any() or set(family.feature_id) != set(recipe['declared_events']):
        raise ValueError('all original declared event definitions required')
    available = {row['event_id'] for row in recipe['events']}
    if not available <= set(family.feature_id) or len(available) != len(recipe['events']):
        raise ValueError('unique declared marker geometries required')
    frames, receipts = [], []
    hashes = {str(path): file_hash(path) for path in (recipe_path, declared_path, definitions_path)}
    seen_barcodes = set()
    for index, packet in enumerate(recipe['bams']):
        if seen_barcodes.intersection(packet['barcode_groups']):
            raise ValueError('disjoint physical-library barcode scopes required')
        seen_barcodes.update(packet['barcode_groups'])
        folder = args.source / f'shard_{index}'
        receipt_path, counts_path = folder / 'manifest.json', folder / 'support.tsv.gz'
        receipt = json.loads(receipt_path.read_text())
        if receipt['complete'] is not True or receipt['shard'] != index or receipt['input'] != packet or receipt['recipe_sha256'] != recipe_hash or receipt['requested_events'] != len(available):
            raise ValueError('all eight complete frozen library count receipts required')
        frame = pd.read_csv(counts_path, sep='\t', dtype={'subject': str})
        if not set(frame.feature_id) <= available or frame.duplicated([*KEYS, 'signature']).any():
            raise ValueError('unique within-library declared signatures required')
        groups = set(map(tuple, packet['barcode_groups'].values()))
        if not set(frame[['subject', 'cell_type', 'primer']].itertuples(index=False, name=None)) <= groups:
            raise ValueError('count rows outside their physical-library cell scope')
        if int(frame['count'].sum()) != receipt['filters']['barcode_umi_event_keys']:
            raise ValueError('molecule-key total differs from complete library receipt')
        hashes.update({str(path): file_hash(path) for path in (receipt_path, counts_path)})
        frames.append(frame)
        receipts.append(receipt)
    if seen_barcodes != set(recipe['barcode_groups']):
        raise ValueError('every retained production cell must have one library')
    signatures = pd.concat(frames, ignore_index=True).groupby([*KEYS, 'signature'], sort=False)['count'].sum().reset_index()
    markers = summarize_markers(signatures)
    genes = gene_key_counts(markers, family.loc[family.feature_id.isin(available)])
    # Deterministic gene partition, never selected by significance or LR labels.
    gene_shards = {gene: index % args.shard_count for index, gene in enumerate(sorted(family.gene_id.unique()))}
    family['shard_index'] = family.gene_id.map(gene_shards)
    markers['shard_index'] = markers.feature_id.map(family.set_index('feature_id').shard_index)
    args.output_dir.mkdir(parents=True)
    family.to_csv(args.output_dir / 'declared_family.tsv.gz', sep='\t', index=False)
    genes.to_csv(args.output_dir / 'gene_keys.tsv.gz', sep='\t', index=False)
    groups = pd.DataFrame(sorted(set(map(tuple, recipe['barcode_groups'].values()))), columns=['subject', 'cell_type', 'primer'])
    groups.to_csv(args.output_dir / 'groups.tsv', sep='\t', index=False)
    shutil.copyfile(definitions_path, args.output_dir / 'feature_definitions.tsv')
    shards = []
    for index in range(args.shard_count):
        folder = args.output_dir / f'shard_{index}'
        folder.mkdir()
        local = markers.loc[markers.shard_index.eq(index)].drop(columns='shard_index')
        path = folder / 'counts.tsv.gz'
        local.to_csv(path, sep='\t', index=False)
        shards.append(dict(shard_index=index, declared_events=int(family.shard_index.eq(index).sum()), count_rows=len(local), counts_sha256=file_hash(path)))
    manifest = dict(input_hashes=hashes, library_receipts=receipts, count_shards=shards, declared_events=len(family), available_marker_geometry=len(available), production_barcodes=len(seen_barcodes), production_cell_qc=recipe['production_cell_qc'], marker_rule='class-unanimous local exon/junction keys, conflicts excluded from either class; separate junction-only ablation', missing_rule='no signature rows means zero observed keys, never removal from declared catalogue', selection=recipe['selection'], scope='complete annotation-catalogue marker counts and deterministic gene partition; no statistical or replication claim', production_changes=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(declared_events=len(family), count_rows=len(markers), genes_with_counts=genes.gene_id.nunique(), shards=len(shards)), indent=2), flush=True)


if __name__ == '__main__':
    main()
