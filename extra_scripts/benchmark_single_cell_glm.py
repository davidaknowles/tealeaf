#!/usr/bin/env python3
"""Benchmark scalable GLM cell-batch sizes on a reproducible cell subset."""

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
    parser.add_argument("--design", required=True, choices=["binary", "weighted"])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--method",
        action="append",
        choices=["factorized", "admm_factorized"],
    )
    parser.add_argument("--batch-cells", action="append", type=int)
    parser.add_argument("--cells", type=int, default=10_000)
    parser.add_argument("--min-cell-umis", type=float, default=500)
    parser.add_argument("--rank", type=int, default=32)
    parser.add_argument("--max-iter", type=int, default=6)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument(
        "--data-backend",
        choices=["auto", "cuda", "pinned"],
        default="auto",
    )
    parser.add_argument(
        "--legacy-api",
        action="store_true",
        help="benchmark the pre-cache solver API from an older PYTHONPATH",
    )
    args = parser.parse_args()
    methods = args.method or ["factorized", "admm_factorized"]
    batch_sizes = args.batch_cells or [8192, 16384, 32768]
    if any(value < 1 for value in batch_sizes):
        raise ValueError("batch sizes must be positive")

    prepared = glm_cv.prepare_alevin_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        ec_design=args.design,
        regularization_target="theta",
    )
    selected = glm_cv.sample_cells_by_count(
        prepared.counts,
        args.cells,
        seed=args.seed,
        min_count=args.min_cell_umis,
        totals=prepared.cell_umi_totals,
    )
    counts = prepared.counts[selected].tocsr()
    torch = glm_solvers._torch()
    rows = []
    for method in methods:
        for batch_cells in batch_sizes:
            if args.device == "cuda":
                torch.cuda.empty_cache()
                torch.cuda.synchronize()
            started = time.perf_counter()
            data = None if args.legacy_api else glm_solvers.SparseGLM(
                counts,
                prepared.compatibility,
                device=args.device,
                batch_cells=batch_cells,
                data_backend=args.data_backend,
            )
            setup_seconds = time.perf_counter() - started
            fit_kwargs = {
                "rank": args.rank,
                "max_iter": args.max_iter,
                "min_iter": args.max_iter,
                "device": args.device,
            }
            if method == "factorized" and not args.legacy_api:
                fit_kwargs.update(minibatch=True, polish_max_iter=1)
            started = time.perf_counter()
            result = (
                glm_solvers.fit_glm(
                    counts, prepared.compatibility, method,
                    batch_cells=batch_cells, **fit_kwargs,
                )
                if args.legacy_api else
                glm_solvers.fit_glm(data, None, method, **fit_kwargs)
            )
            if args.device == "cuda":
                torch.cuda.synchronize()
            fit_seconds = time.perf_counter() - started
            rows.append(
                {
                    "method": method,
                    "batch_cells": int(batch_cells),
                    "data_backend": result.diagnostics.get(
                        "data_backend", "legacy_stream"
                    ),
                    "setup_seconds": setup_seconds,
                    "fit_seconds": fit_seconds,
                    "cells_per_second": result.diagnostics.get(
                        "cells_per_second",
                        len(selected) * result.diagnostics["iterations"]
                        / max(fit_seconds, 1e-12),
                    ),
                    "peak_cuda_memory_bytes": result.diagnostics.get(
                        "peak_cuda_memory_bytes"
                    ),
                    "final_objective_per_cell": (
                        result.diagnostics["objective"][-1]
                        / max(1, int(np.count_nonzero(
                            np.asarray(counts.sum(axis=1)).ravel()
                        )))
                    ),
                    "mean_epoch_seconds": float(np.mean(
                        result.diagnostics.get(
                            "epoch_seconds", [fit_seconds / args.max_iter]
                        )
                    )),
                }
            )
            del result, data
    report = {
        "design": args.design,
        "n_cells": int(len(selected)),
        "min_cell_umis": float(args.min_cell_umis),
        "rank": int(args.rank),
        "iterations": int(args.max_iter),
        "legacy_api": bool(args.legacy_api),
        "results": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
