import numpy as np
import pandas as pd
import pytest

from extra_scripts.audit_ranked_local_marker_scope import diagnose


def fixture():
    ranked = pd.DataFrame(dict(method='m', feature_id=['e', 'unknown', 'e'], contrast_id=['A_B', 'A_B', 'A_C'], rank=[1, 2, 3], p_value=[.001, .01, .1], pooled_replicated=[True, False, False]))
    scope = pd.DataFrame([dict(feature_id='e', marker_variant=variant, event_type='SE', status='ok', outside_marker_sources=2, outside_marker_source_LR_UMIs=5, outside_source_fraction_of_gene=.25, outside_source_fraction_of_marker_source_RNA=.4, event_class_fraction_of_gene=.75) for variant in ('all', 'junction')])
    return ranked, scope


def test_scope_audit_keeps_unknown_events_multiple_contrasts_and_short_prefix():
    ranked, scope = fixture()
    summary, table = diagnose(ranked, scope, 'frozen')
    assert len(table) == 6
    for _, local in table.groupby('marker_variant'):
        assert local.feature_id.tolist() == ranked.feature_id.tolist()
        np.testing.assert_array_equal(local.p_value, ranked.p_value)
        np.testing.assert_array_equal(local['rank'], ranked['rank'])
    assert not summary.complete_prefix.any()
    whole = summary.loc[summary.selection.eq('all original prefix')]
    assert whole.n_original_prefix.eq(3).all() and whole.n_scope_available.eq(2).all() and whole.n_scope_missing.eq(1).all()
    assert whole.events_with_expressed_outside_marker_source.eq(2).all()


def test_scope_audit_never_reranks_by_direction_or_excludes_outside_sources():
    ranked, scope = fixture()
    _, first = diagnose(ranked, scope, 'frozen')
    _, second = diagnose(ranked.assign(pooled_replicated=~ranked.pooled_replicated), scope, 'frozen')
    assert first.feature_id.tolist() == second.feature_id.tolist()
    assert first['rank'].tolist() == second['rank'].tolist()


@pytest.mark.parametrize('problem', ['changed_rank', 'missing_direction', 'duplicate_geometry'])
def test_scope_audit_refuses_changed_prefix_or_ambiguous_scope(problem):
    ranked, scope = fixture()
    if problem == 'changed_rank':
        ranked.loc[0, 'rank'] = 9
    elif problem == 'missing_direction':
        ranked['pooled_replicated'] = ranked.pooled_replicated.astype(object)
        ranked.loc[0, 'pooled_replicated'] = None
    else:
        scope = pd.concat([scope, scope.iloc[:1]])
    with pytest.raises(ValueError):
        diagnose(ranked, scope, 'frozen')
