import pandas as pd

from tealeaf.sc.junction_benchmark import (
    permute_paired_contrasts,
    plan_pairwise_contrasts,
    simes_omnibus,
)


def test_plan_pairwise_contrasts_uses_biological_replicates():
    rows = []
    for subject, condition in [("m1", "A"), ("m2", "A"), ("m3", "B"), ("m4", "B")]:
        for cell_type in ("X", "Y"):
            rows.append((f"{subject}__{cell_type}", cell_type, condition, subject, 10))
    samples = pd.DataFrame(rows, columns=["sample_id", "cell_type", "condition", "subject", "n_cells"])
    contrasts = plan_pairwise_contrasts(samples, min_subjects=2)
    condition = [value for value in contrasts if value["effect"] == "condition"]
    cell_type = [value for value in contrasts if value["effect"] == "cell_type"]
    assert len(condition) == 2
    assert len(cell_type) == 1
    assert cell_type[0]["paired_subjects"] == ["m1", "m2", "m3", "m4"]


def test_paired_contrast_permutation_preserves_subject_pairs():
    samples = pd.DataFrame({
        "sample_id": ["m1_a", "m1_b", "m2_a", "m2_b", "m3_a", "m3_b"],
        "subject": ["m1", "m1", "m2", "m2", "m3", "m3"],
    })
    contrast = {
        "contrast_id": "cell_type__a__b",
        "effect": "cell_type",
        "level_a": "a",
        "level_b": "b",
        "samples_a": ["m1_a", "m2_a", "m3_a"],
        "samples_b": ["m1_b", "m2_b", "m3_b"],
    }
    first = permute_paired_contrasts([contrast], samples, 7)[0]
    second = permute_paired_contrasts([contrast], samples, 7)[0]
    assert first == second
    lookup = samples.set_index("sample_id").subject.to_dict()
    assert [lookup[value] for value in first["samples_a"]] == ["m1", "m2", "m3"]
    assert [lookup[value] for value in first["samples_b"]] == ["m1", "m2", "m3"]
    assert set(first["samples_a"] + first["samples_b"]) == set(samples.sample_id)
    assert first["permutation_seed"] == 7


def test_simes_omnibus_adjusts_features_within_method():
    pairwise = pd.DataFrame(
        {
            "method": ["m", "m", "m"],
            "feature_id": ["a", "a", "b"],
            "p_value": [0.01, 0.04, 0.8],
        }
    )
    result = simes_omnibus(pairwise).set_index("feature_id")
    assert result.loc["a", "p_value"] == 0.02
    assert result.loc["b", "p_value"] == 0.8
