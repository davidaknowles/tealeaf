#!/usr/bin/env python3
"""Benchmark deterministic factorized-GLM polishing from saved factors."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from tealeaf.sc import glm_cv, glm_solvers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument(
        "--design", required=True, choices=("binary", "weighted", "positional")
    )
    parser.add_argument("--initial-prefix", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--inner-steps", type=int, action="append")
    parser.add_argument("--iterations", type=int, default=32)
    parser.add_argument("--batch-cells", type=int, default=16384)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    inner_steps = args.inner_steps or [1, 2, 4, 8]
    prepared = glm_cv.prepare_paired_primer_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        args.primer_pairs,
        ec_design=args.design,
        regularization_target="theta",
        min_eq=5,
        min_half_umis=500,
    )
    factor_path = Path(f"{args.initial_prefix}glm_factors.npz")
    rows = np.loadtxt(f"{args.initial_prefix}glm_rows.txt", dtype=str)
    columns = np.loadtxt(f"{args.initial_prefix}glm_cols.txt", dtype=str)
    if not np.array_equal(rows, prepared.barcodes):
        raise ValueError("initial rows do not match paired preparation")
    if not np.array_equal(columns, prepared.features):
        raise ValueError("initial columns do not match paired preparation")
    with np.load(factor_path) as factors:
        initial = (factors["left"].copy(), factors["right"].copy())

    data = glm_solvers.SparseGLM(
        prepared.counts,
        prepared.compatibility,
        device=args.device,
        batch_cells=args.batch_cells,
        data_backend="auto",
    )
    torch = data.torch
    initial_tensors = glm_solvers._coerce_initial_factors(data, initial)
    initial_objective = data.loss_for_factors(*initial_tensors)
    results = []
    for steps in inner_steps:
        if data.device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(data.device)
            torch.cuda.synchronize(data.device)
        started = time.perf_counter()
        result = glm_solvers.fit_factorized(
            data,
            rank=initial[0].shape[1],
            max_iter=args.iterations,
            min_iter=args.iterations,
            patience=1,
            tol=0.0,
            minibatch=False,
            polish_max_iter=args.iterations,
            exact_inner_steps=steps,
            initial_factors=initial,
            device=args.device,
            batch_cells=args.batch_cells,
        )
        if data.device.type == "cuda":
            torch.cuda.synchronize(data.device)
        elapsed = time.perf_counter() - started
        final_objective = result.diagnostics["objective"][-1]
        results.append({
            "inner_steps": int(steps),
            "iterations": int(args.iterations),
            "elapsed_seconds": elapsed,
            "seconds_per_epoch": elapsed / args.iterations,
            "initial_objective": initial_objective,
            "final_objective": final_objective,
            "objective_reduction": initial_objective - final_objective,
            "reduction_per_second": (
                initial_objective - final_objective
            ) / max(elapsed, 1e-12),
            "last_relative_change": result.diagnostics[
                "objective_relative_change"
            ][-1],
            "peak_cuda_memory_bytes": result.diagnostics[
                "peak_cuda_memory_bytes"
            ],
        })
        print(json.dumps(results[-1], indent=2), flush=True)
        del result
    report = {
        "design": args.design,
        "n_cells": int(prepared.counts.shape[0]),
        "rank": int(initial[0].shape[1]),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2)


if __name__ == "__main__":
    main()
