#!/usr/bin/env python3
"""Fit full-data EC-GLMM nulls or alternatives with selective batching."""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from analyses.benchmark_batched_laplace_lbfgs import (
    filter_metadata_and_counts,
    fit_batched,
    reparameterize_parameters,
)
from extra_scripts.run_ec_block_glmm import (
    local_test_design,
    tensor_data,
    treatment_design,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm, ec_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--batchability", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--model", choices=("null", "alternative"), required=True)
    parser.add_argument("--route", choices=("all", "batched", "scalar"), default="all")
    parser.add_argument(
        "--null-fit-shards",
        help="Glob for completed null-fit directories used to initialize alternatives",
    )
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--minimum-batch-size", type=int, default=4)
    parser.add_argument("--maximum-padding-ratio", type=float, default=2.0)
    parser.add_argument("--maximum-coefficients", type=int, default=128)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--mode-steps", type=int, default=12)
    parser.add_argument("--continuation-mode-steps", type=int, default=30)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--laplace-cache-size", type=int, default=16)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def alternative_fit_id(candidate):
    gene = int(candidate[3])
    rows = np.asarray(candidate[7], dtype=np.int64)
    digest = hashlib.sha1(rows.tobytes()).hexdigest()[:16]
    return f"{gene}|{digest}"


def make_records(candidates, dimensions, model):
    dimension_by_test = dimensions.set_index("test_id").to_dict("index")
    records = []
    seen_alternatives = set()
    for candidate in candidates:
        test_id = candidate[0]
        if test_id not in dimension_by_test:
            continue
        row = dimension_by_test[test_id]
        alternative_id = alternative_fit_id(candidate)
        if model == "alternative":
            if alternative_id in seen_alternatives:
                continue
            seen_alternatives.add(alternative_id)
            fit_id = alternative_id
            coefficients = int(row["alternative_coefficients"])
        else:
            fit_id = test_id
            coefficients = int(row["null_coefficients"])
        records.append(
            {
                "fit_id": fit_id,
                "test_id": test_id,
                "candidate": candidate,
                "isoforms": int(row["isoforms"]),
                "coefficients": coefficients,
                "pair_coefficients": max(
                    int(row["null_coefficients"]),
                    int(row["alternative_coefficients"]),
                ),
                "samples": int(row["samples"]),
                "primer_1_ecs": int(row["primer_1_ecs"]),
                "primer_2_ecs": int(row["primer_2_ecs"]),
                "cost_proxy": int(row["cost_proxy"]),
            }
        )
    return records


def padding_ratio(records):
    actual = sum(
        record["samples"]
        * (record["primer_1_ecs"] + record["primer_2_ecs"])
        * record["isoforms"]
        for record in records
    )
    padded = (
        len(records)
        * max(record["samples"] for record in records)
        * (
            max(record["primer_1_ecs"] for record in records)
            + max(record["primer_2_ecs"] for record in records)
        )
        * records[0]["isoforms"]
    )
    return padded / max(actual, 1)


def build_work_units(records, args):
    grouped = {}
    for record in records:
        grouped.setdefault((record["isoforms"], record["coefficients"]), []).append(
            record
        )
    work = []
    for _, group in grouped.items():
        group.sort(key=lambda record: record["cost_proxy"])
        for start in range(0, len(group), int(args.batch_size)):
            chunk = group[start : start + int(args.batch_size)]
            admitted = (
                len(chunk) >= int(args.minimum_batch_size)
                and max(record["pair_coefficients"] for record in chunk)
                <= int(args.maximum_coefficients)
                and padding_ratio(chunk) <= float(args.maximum_padding_ratio)
            )
            if admitted:
                work.append(
                    {
                        "route": "batched",
                        "records": chunk,
                        "cost": max(record["cost_proxy"] for record in chunk),
                        "padding_ratio": padding_ratio(chunk),
                    }
                )
            else:
                work.extend(
                    {
                        "route": "scalar",
                        "records": [record],
                        "cost": record["cost_proxy"],
                        "padding_ratio": 1.0,
                    }
                    for record in chunk
                )
    return work


def partition_work(work, shard_count):
    shards = [[] for _ in range(int(shard_count))]
    loads = np.zeros(int(shard_count), dtype=float)
    for unit in sorted(work, key=lambda value: value["cost"], reverse=True):
        shard = int(np.argmin(loads))
        shards[shard].append(unit)
        loads[shard] += float(unit["cost"])
    return shards, loads


def prepare_data(record, metadata, counts, gene_ecs, designs, settings, model):
    candidate = record["candidate"]
    _, _, _, gene, transcripts, path_index, _, rows, _, tested_levels = candidate
    local_metadata, nuisance, labels = local_test_design(
        metadata,
        rows,
        tested_levels,
        settings.get("test_effect", "cell_type"),
    )
    local_counts = tuple(matrix[rows] for matrix in counts)
    base, _, _ = local_gene_data(
        local_counts,
        designs,
        transcripts,
        gene_ecs[gene],
        np.ones((len(local_metadata), 1)),
        local_metadata.mouse.to_numpy(),
        drop_zero=False,
    )
    tested = treatment_design(labels, len(tested_levels))
    if model == "alternative":
        tensor = ec_block_glmm.full_fixed_effect_tensor(
            nuisance, tested, len(transcripts) - 1
        )
    else:
        tensor, _, _ = ec_block_glmm.block_fixed_effect_tensors(
            nuisance, tested, path_index
        )
    return tensor_data(base, tensor)


def load_null_initializations(pattern):
    """Load fitted null parameters keyed by the representative test identifier."""
    if pattern is None:
        raise ValueError("alternative fitting requires --null-fit-shards")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise ValueError(f"no null-fit shards matched {pattern}")
    tables = [pd.read_csv(Path(path) / "fits.tsv", sep="\t") for path in paths]
    table = pd.concat(tables, ignore_index=True)
    if table.fit_id.duplicated().any():
        raise ValueError("null-fit shards contain duplicate fit identifiers")
    return {
        str(row.fit_id): np.asarray(json.loads(row.parameters), dtype=float)
        for row in table.itertuples()
    }


def scalar_fit(
    data,
    args,
    initial=None,
    continuation=False,
    objective_cache=None,
):
    mode_steps = (
        int(args.continuation_mode_steps) if continuation else int(args.mode_steps)
    )
    fit = ec_glmm.fit_laplace(
        data,
        initial=initial,
        max_iter=int(args.max_iter),
        mode_steps=mode_steps,
        mode_gradient="implicit",
        objective_cache=objective_cache,
        objective_cache_key=(
            "global_shape",
            ec_glmm.laplace_objective_shape_key(data),
        ),
    )
    attempts = 1
    while not fit["converged"] and attempts <= int(args.retries):
        fit = ec_glmm.fit_laplace(
            data,
            initial=fit["parameters"],
            max_iter=int(args.max_iter),
            mode_steps=int(args.continuation_mode_steps),
            mode_gradient="implicit",
            objective_cache=objective_cache,
            objective_cache_key=(
                "global_shape",
                ec_glmm.laplace_objective_shape_key(data),
            ),
        )
        attempts += 1
    return fit, attempts


def fit_record(
    record,
    fit,
    attempts,
    route,
    unit_seconds,
    padding,
    initialization,
):
    return {
        "fit_id": record["fit_id"],
        "test_id": record["test_id"],
        "route": route,
        "batch_size": 1 if route == "scalar" else np.nan,
        "padding_ratio": padding,
        "converged": bool(fit["converged"]),
        "objective": float(fit["objective"]),
        "iterations": int(fit["iterations"]),
        "scaled_gradient_norm": float(fit["scaled_gradient_norm"]),
        "mode_score_norm": float(fit["mode_score_norm"]),
        "attempts": int(attempts),
        "unit_seconds": float(unit_seconds),
        "parameters": json.dumps(np.asarray(fit["parameters"], dtype=float).tolist()),
        "initialization": initialization,
    }


def run_unit(
    unit,
    metadata,
    counts,
    gene_ecs,
    designs,
    settings,
    args,
    null_initializations,
    objective_cache,
):
    records = unit["records"]
    data = [
        prepare_data(
            record,
            metadata,
            counts,
            gene_ecs,
            designs,
            settings,
            args.model,
        )
        for record in records
    ]
    initial = None
    if args.model == "alternative":
        initial = []
        for record, alternative_data in zip(records, data):
            null_data = prepare_data(
                record,
                metadata,
                counts,
                gene_ecs,
                designs,
                settings,
                "null",
            )
            if record["test_id"] not in null_initializations:
                raise ValueError(
                    f"missing fitted null for alternative {record['fit_id']}"
                )
            initial.append(
                reparameterize_parameters(
                    null_initializations[record["test_id"]],
                    null_data.fixed_effect_tensor,
                    alternative_data.fixed_effect_tensor,
                )
            )
    started = time.perf_counter()
    initialization = (
        "fitted_null_projection" if initial is not None else "default"
    )
    if unit["route"] == "scalar":
        fit, attempts = scalar_fit(
            data[0],
            args,
            initial=None if initial is None else initial[0],
            objective_cache=objective_cache,
        )
        elapsed = time.perf_counter() - started
        return [
            fit_record(
                records[0],
                fit,
                attempts,
                "scalar",
                elapsed,
                1.0,
                initialization,
            )
        ]
    result, _, gradient_norm, mode_score, converged = fit_batched(
        data,
        int(args.max_iter),
        int(args.mode_steps),
        initial=initial,
    )
    elapsed = time.perf_counter() - started
    rows = []
    for index, record in enumerate(records):
        if converged[index]:
            rows.append(
                {
                    "fit_id": record["fit_id"],
                    "test_id": record["test_id"],
                    "route": "batched",
                    "batch_size": len(records),
                    "padding_ratio": unit["padding_ratio"],
                    "converged": True,
                    "objective": float(result.values[index]),
                    "iterations": int(result.iterations[index]),
                    "scaled_gradient_norm": float(gradient_norm[index]),
                    "mode_score_norm": float(mode_score[index]),
                    "attempts": 1,
                    "unit_seconds": elapsed,
                    "parameters": json.dumps(
                        np.asarray(result.parameters[index], dtype=float).tolist()
                    ),
                    "initialization": initialization,
                }
            )
            continue
        fit, attempts = scalar_fit(
            data[index],
            args,
            initial=result.parameters[index],
            continuation=True,
            objective_cache=objective_cache,
        )
        route = "batched_scalar_continuation"
        record_initialization = initialization
        if not fit["converged"]:
            rescue, rescue_attempts = scalar_fit(
                data[index],
                args,
                initial=None if initial is None else initial[index],
                objective_cache=objective_cache,
            )
            attempts += rescue_attempts
            route = "batched_scalar_rescue"
            record_initialization = f"{initialization}_rescue"
            if rescue["converged"] or rescue["objective"] < fit["objective"]:
                fit = rescue
        rows.append(
            fit_record(
                record,
                fit,
                attempts,
                route,
                elapsed,
                unit["padding_ratio"],
                record_initialization,
            )
        )
    return rows


def main():
    args = parse_args()
    with args.data_cache.open("rb") as handle:
        groups, counts, _, _, gene_ecs, designs = pickle.load(handle)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached.get("settings", {})
    metadata, counts = filter_metadata_and_counts(groups, counts, settings)
    dimensions = pd.read_csv(args.batchability, sep="\t")
    records = make_records(cached["candidates"], dimensions, args.model)
    work = build_work_units(records, args)
    if args.route != "all":
        work = [unit for unit in work if unit["route"] == args.route]
    shards, loads = partition_work(work, args.shard_count)
    assigned = shards[int(args.shard_index)]
    if args.dry_run:
        summary = {
            "model": args.model,
            "fits": len(records),
            "work_units": len(work),
            "batched_fits": sum(
                len(unit["records"]) for unit in work if unit["route"] == "batched"
            ),
            "scalar_fits": sum(
                len(unit["records"]) for unit in work if unit["route"] == "scalar"
            ),
            "shard_load_minimum": float(loads.min()),
            "shard_load_maximum": float(loads.max()),
        }
        print(json.dumps(summary, indent=2))
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    null_initializations = (
        load_null_initializations(args.null_fit_shards)
        if args.model == "alternative"
        else None
    )
    started = time.perf_counter()
    results = []
    objective_cache = {}
    for unit in assigned:
        results.extend(
            run_unit(
                unit,
                metadata,
                counts,
                gene_ecs,
                designs,
                settings,
                args,
                null_initializations,
                objective_cache,
            )
        )
        if unit["route"] == "batched":
            try:
                import jax

                jax.clear_caches()
            except ImportError:
                pass
        while len(objective_cache) > int(args.laplace_cache_size):
            objective_cache.pop(next(iter(objective_cache)))
    pd.DataFrame(results).to_csv(args.output_dir / "fits.tsv", sep="\t", index=False)
    summary = {
        "model": args.model,
        "route": args.route,
        "shard_index": int(args.shard_index),
        "shard_count": int(args.shard_count),
        "work_units": len(assigned),
        "fits": len(results),
        "seconds": time.perf_counter() - started,
        "routes": pd.Series([result["route"] for result in results])
        .value_counts()
        .to_dict(),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
