import json
import pickle
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from extra_scripts.prepare_sequence_event_control import main


def fixture(tmp_path):
    source = tmp_path / "inputs/original_binary_paired"
    source.mkdir(parents=True)
    designs = (csr_matrix([[1., 0], [2., 3.]]), csr_matrix([[4., 0], [5., 6.]]))
    counts = (csr_matrix([[3, 7]]), csr_matrix([[11, 13]]))
    prepared = (["subject"], counts, ["gene"], [np.array([0, 1])], [np.array([0, 1])], designs)
    with (source / "prepared.pkl").open("wb") as handle:
        pickle.dump(prepared, handle)
    (source / "features.txt").write_text("tx1\ntx2\n")
    shard = tmp_path / "kernels/original_binary/shard_0"
    shard.mkdir(parents=True)
    settings = dict(source="original_binary", shard_count=1, read_lengths=[150], terminal_start_windows=["full"], background_fraction=.01, selection="whole screened gene catalog", requested_genes=["gene"])
    (shard / "settings.json").write_text(json.dumps(settings))
    pd.DataFrame([dict(gene_id="gene", primer=primer, source="original_binary", status="ok", full_positive_support=True, retained_molecule_fraction=1.) for primer in (0, 1)]).to_csv(shard / "genes.tsv.gz", sep="\t", index=False)
    dt = np.array([[.2, 0], [.8, 1.]])
    packet = shard / "gene_L150_Wfull_E0.01.npz"
    np.savez(packet, ec_indices=[0, 1], transcripts=["tx1", "tx2"], background_fraction=.01, dt=dt, rh=dt * [10, 20])
    return prepared, packet


def invoke(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["prepare", "--input-root", str(tmp_path / "inputs"), "--kernel-root", str(tmp_path / "kernels"), "--source", "original_binary", "--output-dir", str(tmp_path / "output"), "--shard-count", "1"])
    main()


def test_preparation_preserves_counts_order_and_support(tmp_path, monkeypatch):
    before, _ = fixture(tmp_path)
    invoke(tmp_path, monkeypatch)
    with (tmp_path / "output/prepared.pkl").open("rb") as handle:
        after = pickle.load(handle)
    for old, new in zip(before[1], after[1]):
        np.testing.assert_array_equal(old.toarray(), new.toarray())
    np.testing.assert_array_equal(after[5][0].toarray(), [[.2, 0], [.8, 1]])
    for old, new in zip(before[5], after[5]):
        np.testing.assert_array_equal(old.indices, new.indices)
        np.testing.assert_array_equal(old.indptr, new.indptr)
    manifest = json.loads((tmp_path / "output/manifest.json").read_text())
    assert manifest["patched_genes"] == 1 and not manifest["production_changes"]
    with pytest.raises(FileExistsError):
        invoke(tmp_path, monkeypatch)


@pytest.mark.parametrize("defect", ["alignment", "support"])
def test_invalid_packets_do_not_produce_prepared_input(tmp_path, monkeypatch, defect):
    _, packet = fixture(tmp_path)
    with np.load(packet) as loaded:
        values = dict(loaded)
    if defect == "alignment":
        values["transcripts"] = np.array(["tx2", "tx1"])
    else:
        values["dt"][0, 1] = .1
    np.savez(packet, **values)
    with pytest.raises(ValueError):
        invoke(tmp_path, monkeypatch)
    assert not (tmp_path / "output/prepared.pkl").exists()
