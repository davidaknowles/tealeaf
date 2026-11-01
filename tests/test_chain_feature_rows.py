import json
from types import SimpleNamespace

import pandas as pd

from extra_scripts.assess_paired_inference_audit import chain_feature_rows


def test_chain_feature_rows_uses_fit_path_members():
    # Signatures 0 and 1 have no junction and form one group in the fit; regrouping the block paths (where path 1 has a
    # junction) would split them and shift every later tested path by one.
    signatures = [[[100, 150]], [[100, 170]], [[300, 400]], [[500, 600]]]
    paths = [[[0, 150]], [[0, 50], [80, 170]], [[0, 50], [300, 400]], [[0, 50], [500, 600]]]
    block = {"signatures": signatures, "paths": paths}
    row = SimpleNamespace(kept_paths=json.dumps([1, 2]), path_members=json.dumps([[signatures[2]], [signatures[3]]]))
    features = pd.DataFrame({"transcript_id": ["a", "b"]})
    chains = {"a": [(10, 50), (300, 390)], "b": [(10, 50), (500, 590)]}
    local = chain_feature_rows(block, row, features, chains)
    assert dict(zip(local.transcript_id, local.path_number)) == {"a": 1, "b": 2}
