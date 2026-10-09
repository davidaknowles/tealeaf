import pandas as pd
import pytest

from extra_scripts.audit_lr_catalog_retention import catalog_retention


def test_catalog_audit_retains_every_prefix_identity_without_reranking():
    ranked = pd.DataFrame(dict(feature_id=["SUPPA2:g1.2;old", "SUPPA2:g1.2;new", "SUPPA2:g2.1;new"], rank=[1, 2, 3], contrast_id=["A__B"] * 3))
    reference = pd.DataFrame(dict(feature_id=["SUPPA2:g1.2;old"], gene_id=["g1"]))
    definitions = pd.DataFrame(dict(feature_id=ranked.feature_id, status=["ok", "ok", "unavailable"]))
    original = ranked.copy(deep=True)
    output = catalog_retention(ranked, reference, definitions)
    assert output.feature_id.tolist() == ranked.feature_id.tolist()
    assert output['rank'].tolist() == ranked['rank'].tolist()
    assert output.prior_definition_supported.tolist() == [True, False, False]
    assert output.prior_gene_supported.tolist() == [True, True, False]
    assert output.full_catalog_geometry_status.tolist() == ["ok", "ok", "unavailable"]
    pd.testing.assert_frame_equal(ranked, original)


def test_catalog_audit_rejects_ambiguous_reference_definitions():
    ranked = pd.DataFrame(dict(feature_id=["SUPPA2:g1;event"], rank=[1]))
    reference = pd.DataFrame(dict(feature_id=["e", "e"], gene_id=["g1", "g1"]))
    definitions = pd.DataFrame(dict(feature_id=["e"], status=["ok"]))
    with pytest.raises(ValueError):
        catalog_retention(ranked, reference, definitions)
