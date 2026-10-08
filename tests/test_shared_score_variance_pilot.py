import numpy as np
import pandas as pd

from extra_scripts.audit_shared_score_variance import freeze_training_cases, prepare_panel, signed_values


def test_training_selection_precedes_fit_availability_and_is_order_invariant():
    table = pd.DataFrame(dict(test_id=["a", "b", "c", "d"], gene_id=["g1", "g1", "g2", "g2"], p_value=[.01, .7, .0001, .5], complete_subject_fits=[True] * 4))
    frozen = freeze_training_cases(table)
    revised = table.copy()
    revised["p_value"] = revised.p_value[::-1].to_numpy()
    revised.loc[revised.test_id.isin(frozen.test_id), "complete_subject_fits"] = False
    selected = freeze_training_cases(revised.sample(frac=1., random_state=57))
    assert selected.test_id.tolist() == frozen.test_id.tolist()
    assert not selected.complete_subject_fits.any()
    assert len(selected) == 2


def test_signed_values_use_original_subject_rng_before_rank_masking():
    import zlib
    ids = ["a", "b"]
    records = {key: pd.DataFrame(dict(subject=[f"s{x}" for x in range(6)], score=np.arange(6) + 1., information=[1., 1., 1e-14, 1., 1., 1.], biological_shape=[4.] * 6, reference_information=[1.] * 6)) for key in ids}
    panel, table = prepare_panel(records, ids)
    for replicate in range(32):
        source = np.concatenate([np.random.default_rng(np.random.SeedSequence((20260927, zlib.crc32(key.encode()), replicate))).choice((-1., 1.), size=6) for key in ids])
        np.testing.assert_array_equal(signed_values(panel, table, ids, 20260927, replicate), panel.values * source[panel.source_positions])
    assert panel.n_subjects.tolist() == [5, 5]
