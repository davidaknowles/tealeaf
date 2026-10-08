import numpy as np
import pandas as pd
import pytest

from extra_scripts.pilot_library_read_models import count_tensor, marker_lookup
from extra_scripts.assess_library_read_model_pilot import validate_family


def test_marker_ablation_retains_conflicts_without_treating_exons_as_junctions():
    frame = pd.DataFrame(dict(feature_id=["e"] * 5, subject=["s"] * 5, cell_type=["B"] * 5, primer=["poly(dT)"] * 5, signature=[17, 8, 25, 12, 4], count=[2, 3, 5, 7, 11]))
    original = frame.copy(deep=True)
    all_keys = marker_lookup(frame, junction_only=False)
    junctions = marker_lookup(frame, junction_only=True)
    assert all_keys[("e", "s", "B", "poly(dT)")] == dict(included=13, excluded=3)
    assert junctions[("e", "s", "B", "poly(dT)")] == dict(included=11, excluded=8)
    values = count_tensor("e", ["s", "unmeasured"], ["A", "B"], junctions)
    assert values.shape == (2, 2, 2, 2)
    np.testing.assert_array_equal(values[0, 0, 1], [11, 8])
    assert not values[0, :, 0].any()
    assert not values[1].any()
    assert not values[:, 1].any()
    pd.testing.assert_frame_equal(frame, original)


@pytest.mark.parametrize("subjects,levels", [(["s", "s"], ["A", "B"]), (["s"], ["A", "A"]), (["s"], ["A"])])
def test_pilot_tensor_rejects_duplicated_subjects_or_wrong_contrast(subjects, levels):
    with pytest.raises(ValueError):
        count_tensor("e", subjects, levels, {})


def fitted_family():
    cases = pd.DataFrame([dict(fold=0, test_id="test")])
    table = pd.DataFrame([dict(fold=0, test_id="test", model=model, variant=variant, counts_sha256=variant, p_value=1., converged=False) for model in ("conditional", "unconditional") for variant in ("all local markers", "junction markers")])
    return table, cases


def test_pilot_assessor_keeps_every_failed_model_variant_identity():
    table, cases = fitted_family()
    validate_family(table, cases)


@pytest.mark.parametrize("problem", ["missing", "duplicate", "changed_counts", "failed_p"])
def test_pilot_assessor_rejects_incomplete_family_or_nonmatched_inputs(problem):
    table, cases = fitted_family()
    if problem == "missing":
        table = table.iloc[:-1]
    elif problem == "duplicate":
        table = pd.concat([table, table.iloc[:1]])
    elif problem == "changed_counts":
        table.loc[0, "counts_sha256"] = "changed"
    else:
        table.loc[0, "p_value"] = .01
    with pytest.raises(ValueError):
        validate_family(table, cases)
