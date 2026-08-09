"""Tests for the selective full-data batching schedule."""

from types import SimpleNamespace

from analyses.run_hybrid_batched_ec_glmm import (
    build_work_units,
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
