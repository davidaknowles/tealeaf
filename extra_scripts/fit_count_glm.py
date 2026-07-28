#!/usr/bin/env python3
"""Fit the count-aware flat transcript GLM to paired-primer data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from tealeaf.sc import count_glm, glm_cv, hierarchical


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument("--transcript-to-gene", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument(
        "--ec-design",
        choices=["binary", "weighted", "positional"],
        default="positional",
    )
    parser.add_argument(
        "--primer-sampling-model",
        choices=["effective_length", "oligodt_tpm", "all_tpm"],
        default="oligodt_tpm",
    )
    parser.add_argument("--rank", type=int, default=32)
    parser.add_argument("--hvg", type=int, default=2_000)
    parser.add_argument("--min-feature-cells", type=int, default=200)
    parser.add_argument(
        "--gene-family",
        choices=["negative_binomial", "poisson", "standardized_log"],
        default="negative_binomial",
    )
    parser.add_argument("--concentration", type=float, default=10.0)
    parser.add_argument("--initialization", choices=["pca", "random"], default="pca")
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-cell-umis", type=float, default=500)
    parser.add_argument("--min-eq", type=int, default=5)
    parser.add_argument("--max-cells", type=int)
    parser.add_argument("--gene-epochs", type=int, default=10)
    parser.add_argument("--ec-epochs", type=int, default=10)
    parser.add_argument("--joint-epochs", type=int, default=10)
    parser.add_argument("--batch-cells", type=int, default=256)
    parser.add_argument("--gene-learning-rate", type=float, default=1e-2)
    parser.add_argument("--ec-learning-rate", type=float, default=3e-3)
    parser.add_argument("--joint-learning-rate", type=float, default=1e-3)
    parser.add_argument("--gene-weight", type=float, default=1.0)
    parser.add_argument("--loading-penalty", type=float, default=1e-4)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    started = time.perf_counter()
    print(json.dumps({"event": "paired_preparation_started"}), flush=True)
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
    print(
        json.dumps(
            {
                "event": "paired_preparation_complete",
                "seconds": time.perf_counter() - started,
            }
        ),
        flush=True,
    )
    selected_cells = glm_cv.sample_cells_by_count(
        prepared.counts,
        0,
        min_count=args.min_cell_umis,
        totals=prepared.cell_umi_totals,
    )
    if args.max_cells is not None and len(selected_cells) > args.max_cells:
        rng = np.random.default_rng(args.seed)
        selected_cells = np.sort(
            rng.choice(selected_cells, size=args.max_cells, replace=False)
        )
    projection_started = time.perf_counter()
    print(
        json.dumps(
            {
                "event": "gene_projection_started",
                "selected_cells": int(len(selected_cells)),
            }
        ),
        flush=True,
    )
    data = hierarchical.prepare_hierarchical_data(
        prepared, selected_cells, args.transcript_to_gene
    )
    print(
        json.dumps(
            {
                "event": "data_prepared",
                "seconds": time.perf_counter() - started,
                "gene_projection_seconds": (time.perf_counter() - projection_started),
                "cells": len(data.barcodes),
                "genes": len(data.genes),
                "transcripts": len(data.transcripts),
            }
        ),
        flush=True,
    )

    initial_state, selected_genes, _, pca_diagnostics = (
        hierarchical.log_gene_randomized_pca(
            data.gene_counts,
            args.rank,
            n_hvg=args.hvg,
            min_feature_cells=args.min_feature_cells,
            seed=args.seed,
        )
    )
    print(
        json.dumps({"event": "gene_pca_complete", **pca_diagnostics}),
        flush=True,
    )
    model = count_glm.CountAwareGLM(
        data,
        args.rank,
        selected_genes,
        gene_family=args.gene_family,
        concentration=args.concentration,
        device=args.device,
        initial_state=(initial_state if args.initialization == "pca" else None),
        seed=args.seed,
    )
    model.write(
        f"{args.output_prefix}{args.initialization}_",
        diagnostics={
            "stage": args.initialization,
            "feature_selection": "log_gene_variance",
            **pca_diagnostics,
        },
    )

    histories = {}
    specs = (
        ("gene", args.gene_epochs, args.gene_learning_rate),
        ("ec", args.ec_epochs, args.ec_learning_rate),
        ("joint", args.joint_epochs, args.joint_learning_rate),
    )
    for stage, epochs, learning_rate in specs:
        if epochs <= 0:
            continue
        histories[stage] = model.fit_stage(
            stage,
            epochs=epochs,
            batch_cells=args.batch_cells,
            learning_rate=learning_rate,
            gene_weight=args.gene_weight,
            loading_penalty=args.loading_penalty,
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
                "gene_weight": args.gene_weight,
                "loading_penalty": args.loading_penalty,
                **pca_diagnostics,
            },
        )

    Path(f"{args.output_prefix}run_summary.json").write_text(
        json.dumps(
            {
                "cells": int(len(data.barcodes)),
                "rank": args.rank,
                "initialization": args.initialization,
                "gene_family": args.gene_family,
                "concentration": args.concentration,
                "ec_design": args.ec_design,
                "primer_sampling_model": args.primer_sampling_model,
                "selected_genes": int(len(selected_genes)),
                "histories": histories,
                "pca": pca_diagnostics,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
