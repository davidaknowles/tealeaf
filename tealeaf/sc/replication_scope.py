"""Scope diagnostics for fixed ranked associations, never a new endpoint."""

import json

import numpy as np
import pandas as pd

from .replication_audit import ranked_direction_summary


def diagnose_ranked_marker_scope(ranked, scope, label):
    """Keep original ranks/unknown events, summarize scope without reranking."""
    ranked_direction_summary(ranked)
    columns = ['feature_id', 'event_type', 'status', 'outside_marker_sources', 'outside_marker_source_LR_UMIs', 'outside_source_fraction_of_gene', 'outside_source_fraction_of_marker_source_RNA', 'event_class_fraction_of_gene']
    if scope.duplicated(['feature_id', 'marker_variant']).any():
        raise ValueError('unique original event/marker scope required')
    records, frames = [], []
    for variant, local_scope in scope.groupby('marker_variant'):
        table = ranked.merge(local_scope[columns], on='feature_id', how='left', validate='many_to_one', suffixes=('', '_scope'), indicator=True)
        table['scope_available'] = table._merge.eq('both') & table.status.eq('ok')
        table['source_label'] = label
        table['marker_variant'] = variant
        table['scope_event_type'] = table.event_type_scope if 'event_type_scope' in table else table.event_type
        frames.append(table.drop(columns='_merge'))
        for method, group in table.groupby('method'):
            for cutoff in (100, 200):
                prefix = group.loc[group['rank'].le(cutoff)]
                selections = (('all original prefix', prefix), ('agreeing original prefix', prefix.loc[prefix.pooled_replicated.astype(str).str.lower().eq('true')]), ('disagreeing original prefix', prefix.loc[prefix.pooled_replicated.astype(str).str.lower().eq('false')]))
                for selection, selected in selections:
                    available = selected.loc[selected.scope_available]
                    records.append(dict(source_label=label, method=method, marker_variant=variant, cutoff=cutoff, complete_prefix=len(prefix) == cutoff, selection=selection, n_original_prefix=len(prefix), n_selected=len(selected), n_scope_available=len(available), n_scope_missing=len(selected) - len(available), events_with_outside_marker_source=int(available.outside_marker_sources.gt(0).sum()), events_with_expressed_outside_marker_source=int(available.outside_marker_source_LR_UMIs.gt(0).sum()), median_outside_source_fraction_of_gene=float(available.outside_source_fraction_of_gene.median()) if len(available) else np.nan, median_outside_source_fraction_of_marker_source_RNA=float(available.outside_source_fraction_of_marker_source_RNA.median()) if len(available) else np.nan, median_event_class_fraction_of_gene=float(available.event_class_fraction_of_gene.median()) if len(available) else np.nan, event_types=json.dumps(available.scope_event_type.value_counts().to_dict(), sort_keys=True)))
    return pd.DataFrame(records), pd.concat(frames, ignore_index=True)
