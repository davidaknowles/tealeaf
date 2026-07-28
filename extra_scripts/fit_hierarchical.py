#!/usr/bin/env python3
"""Fit the paired-primer hierarchical gene/isoform multinomial model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from tealeaf.sc import glm_cv, hierarchical


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument("--transcript-to-gene", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--ec-design", default="positional")
    parser.add_argument("--primer-sampling-model", default="oligodt_tpm")
    parser.add_argument("--min-eq", type=int, default=5)
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-cell-umis", type=float, default=500)
    parser.add_argument(
        "--max-cells",
        type=int,
        help="optional deterministic cell cap for smoke tests",
    )
    parser.add_argument("--rank", type=int, default=32)
    parser.add_argument("--hvg", type=int, default=2_000)
    parser.add_argument("--pca-power-iterations", type=int, default=2)
    parser.add_argument("--gene-epochs", type=int, default=10)
    parser.add_argument("--isoform-epochs", type=int, default=10)
    parser.add_argument("--joint-epochs", type=int, default=10)
    parser.add_argument("--batch-cells", type=int, default=256)
    parser.add_argument("--gene-learning-rate", type=float, default=1e-2)
    parser.add_argument("--isoform-learning-rate", type=float, default=3e-3)
    parser.add_argument("--joint-learning-rate", type=float, default=1e-3)
    parser.add_argument("--isoform-penalty", type=float, default=1e-4)
    parser.add_argument(
        "--joint-gene-weight",
        type=float,
        default=1.0,
        help="weight on gene multinomial loss during joint fitting",
    )
    parser.add_argument("--device", default="auto")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--write-pca-initialization",
        action="store_true",
        help="write OUTPUT_PREFIX gene_pca_* files before optimization",
    )
    args = parser.parse_args()

    started = time.perf_counter()
    prepared = glm_cv.prepare_paired_primer_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        args.primer_pairs,
        ec_design=args.ec_design,
        regularization_target="theta",
        min_eq=args.min_eq,
        min_half_umis=args.min_half_umis,
        primer_sampling_model=args.primer_sampling_model,
    )
    selected = glm_cv.sample_cells_by_count(
        prepared.counts,
        0,
        min_count=args.min_cell_umis,
        totals=prepared.cell_umi_totals,
    )
    if args.max_cells is not None and len(selected) > args.max_cells:
        rng = np.random.default_rng(args.seed)
        selected = np.sort(rng.choice(selected, size=args.max_cells, replace=False))
    data = hierarchical.prepare_hierarchical_data(
        prepared, selected, args.transcript_to_gene
    )
    print(
        json.dumps(
            {
                "event": "data_prepared",
                "seconds": time.perf_counter() - started,
                "cells": len(data.barcodes),
                "genes": len(data.genes),
                "transcripts": len(data.transcripts),
                "gene_count_nonzeros": int(data.gene_counts.nnz),
                "primer_count_nonzeros": [int(matrix.nnz) for matrix in data.counts],
            }
        ),
        flush=True,
    )

    pca_started = time.perf_counter()
    initial_state, _, _, pca_diagnostics = hierarchical.log_gene_randomized_pca(
        data.gene_counts,
        args.rank,
        n_hvg=args.hvg,
        power_iterations=args.pca_power_iterations,
        seed=args.seed,
    )
    pca_diagnostics["gene_pca_seconds"] = time.perf_counter() - pca_started
    print(
        json.dumps({"event": "gene_pca_complete", **pca_diagnostics}),
        flush=True,
    )
    model = hierarchical.HierarchicalModel(
        data,
        args.rank,
        device=args.device,
        initial_state=initial_state,
        seed=args.seed,
    )
    if args.write_pca_initialization:
        model.write(
            f"{args.output_prefix}gene_pca_",
            diagnostics={"stage": "gene_pca", **pca_diagnostics},
        )

    histories = {}
    stage_specs = (
        ("gene", args.gene_epochs, args.gene_learning_rate),
        ("isoform", args.isoform_epochs, args.isoform_learning_rate),
        ("joint", args.joint_epochs, args.joint_learning_rate),
    )
    for stage, epochs, learning_rate in stage_specs:
        if epochs <= 0:
            continue
        histories[stage] = model.fit_stage(
            stage,
            epochs=epochs,
            batch_cells=args.batch_cells,
            learning_rate=learning_rate,
            isoform_penalty=args.isoform_penalty,
            joint_gene_weight=args.joint_gene_weight,
            seed=args.seed,
            progress_callback=lambda row: print(
                json.dumps({"event": "epoch_complete", **row}), flush=True
            ),
        )
        model.write(
            f"{args.output_prefix}{stage}_",
            diagnostics={
                "stage": stage,
                "histories": histories,
                "isoform_penalty": args.isoform_penalty,
                "joint_gene_weight": args.joint_gene_weight,
                **pca_diagnostics,
            },
        )

    final_stage = next(reversed(histories), "gene_pca")
    Path(f"{args.output_prefix}run_summary.json").write_text(
        json.dumps(
            {
                "final_stage": final_stage,
                "selected_cells": int(len(selected)),
                "rank": args.rank,
                "histories": histories,
                "pca": pca_diagnostics,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
