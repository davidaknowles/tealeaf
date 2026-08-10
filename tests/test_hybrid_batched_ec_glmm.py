"""Tests for the selective full-data batching schedule."""

from types import SimpleNamespace
import json

import numpy as np
import pandas as pd

from analyses.merge_hybrid_batched_ec_glmm import bh_over_converged
from analyses.run_hybrid_batched_ec_glmm import (
    build_work_units,
    load_null_initializations,
    partition_work,
)


def record(index, *, samples=20, ecs=10, coefficients=12, pair_coefficients=20):
    return {
        "fit_id": str(index),
        "test_id": str(index),
        "candidate": None,
        "isoforms": 3,
        "coefficients": coefficients,
        "pair_coefficients": pair_coefficients,
        "samples": samples,
        "primer_1_ecs": ecs,
        "primer_2_ecs": ecs,
        "cost_proxy": samples * ecs,
    }


def settings():
    return SimpleNamespace(
        batch_size=8,
        minimum_batch_size=4,
        maximum_padding_ratio=2.0,
        maximum_coefficients=128,
    )


def test_scheduler_admits_moderate_low_padding_batch():
    work = build_work_units([record(index) for index in range(8)], settings())
    assert len(work) == 1
    assert work[0]["route"] == "batched"
    assert len(work[0]["records"]) == 8


def test_scheduler_rejects_large_pair_and_high_padding():
    large = build_work_units(
        [record(index, pair_coefficients=129) for index in range(8)],
        settings(),
    )
    assert len(large) == 8
    assert {unit["route"] for unit in large} == {"scalar"}
    padded_records = [record(index, samples=5) for index in range(7)]
    padded_records.append(record(7, samples=100))
    padded = build_work_units(padded_records, settings())
    assert len(padded) == 8
    assert {unit["route"] for unit in padded} == {"scalar"}


def test_work_partition_preserves_every_unit():
    work = build_work_units([record(index) for index in range(24)], settings())
    shards, loads = partition_work(work, 2)
    assert sum(len(shard) for shard in shards) == len(work)
    assert (loads > 0).all()


def test_bh_over_converged_excludes_failed_fits():
    adjusted = bh_over_converged(
        [0.01, 0.02, 0.03, np.nan],
        [True, True, False, False],
    )
    np.testing.assert_allclose(adjusted, [0.02, 0.02, 1.0, 1.0])


def test_load_null_initializations_reads_parameter_vectors(tmp_path):
    first = tmp_path / "null_scalar_shard0"
    second = tmp_path / "null_batched_shard0"
    first.mkdir()
    second.mkdir()
    pd.DataFrame(
        {"fit_id": ["a"], "parameters": [json.dumps([1.0, 2.0])]}
    ).to_csv(first / "fits.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"fit_id": ["b"], "parameters": [json.dumps([3.0, 4.0])]}
    ).to_csv(second / "fits.tsv", sep="\t", index=False)
    loaded = load_null_initializations(str(tmp_path / "null_*_shard*"))
    np.testing.assert_array_equal(loaded["a"], [1.0, 2.0])
    np.testing.assert_array_equal(loaded["b"], [3.0, 4.0])
