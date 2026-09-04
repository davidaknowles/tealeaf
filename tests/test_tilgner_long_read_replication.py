import numpy as np
import pandas as pd

from extra_scripts.assess_tilgner_long_read_replication import biological_replicate, bh_adjust, normalized_difference, summarize_replication, tilgner_cell_type, vector_agreement, wilson_interval


def test_tilgner_cell_type_mapping():
    assert tilgner_cell_type("P56", "Hippocampus", "ExciteNeuron", "ExciteDG") == "EX_hippocampus_granule"
    assert tilgner_cell_type("P56", "Striatum", "InhibNeuron", "D1MSN") == "INH_medium_spiny"
    assert tilgner_cell_type("P28", "Hippocampus", "ExciteNeuron", "ExciteDG") is None
    assert tilgner_cell_type("P56", "Hippocampus", "Oligo", "OPCs") is None
    assert tilgner_cell_type("P56", "Hippocampus", "Oligo", "OPCs", "all") == "ODC"


def test_replicate_mapping():
    assert biological_replicate("Hippocampus", "M1") == 1
    assert biological_replicate("VisCortex", "M9") == 2
    assert biological_replicate("Hippocampus", "M8") is None


def test_vector_direction_and_agreement():
    delta = normalized_difference(np.array([8, 2]), np.array([2, 8]))
    assert np.allclose(delta, [-0.6, 0.6])
    dot, cosine = vector_agreement(np.array([-0.4, 0.4]), delta)
    assert dot > 0
    assert np.isclose(cosine, 1)


def test_bh_and_wilson_helpers():
    assert np.allclose(bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
    low, high = wilson_interval(7, 10)
    assert low < 0.7 < high


def test_strict_replication_uses_conditional_sign_flip_null():
    results = pd.DataFrame({"mapping_complete": [True] * 8, "minimum_pooled_depth": [20] * 8, "minimum_replicate_depth": [10] * 8, "original_effect_norm": [0.1] * 8, "pooled_replicated": [True] * 5 + [False] * 3, "both_replicates_replicated": [True] * 3 + [False] * 5, "replicate_1_dot_product": [1, 1, 1, -1, -1, -1, 1, -1], "replicate_2_dot_product": [1, 1, 1, -1, -1, 1, -1, 1], "technical_g_test_p_value": [1.0] * 8})
    summary = summarize_replication(results)
    selected = summary[(summary["minimum_depth"] == 20) & (summary["minimum_original_effect_norm"] == 0.1)]
    assert dict(zip(selected["endpoint"], selected["null_rate"])) == {"pooled direction": 0.5, "both biological replicates": 5 / 16}
