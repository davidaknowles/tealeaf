#!/usr/bin/env python3
"""Tune a scalable single-cell GLM hyperparameter by held-out molecule counts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from tealeaf.sc import glm_cv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument(
        "--design", required=True, choices=["binary", "weighted", "positional"]
    )
    parser.add_argument(
        "--method",
        required=True,
        choices=["factorized", "admm_factorized", "frank_wolfe_penalized"],
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument(
        "--primer-sampling-model",
        choices=["effective_length", "oligodt_tpm", "all_tpm"],
        default="effective_length",
    )
    parser.add_argument("--multiplier", action="append", type=float)
    parser.add_argument("--rank-candidate", action="append", type=int)
    parser.add_argument("--max-rank", type=int, default=256)
    parser.add_argument(
        "--cells",
        type=int,
        default=5_000,
        help="number of eligible CV cells to sample; use 0 for all",
    )
    parser.add_argument(
        "--min-cell-umis",
        type=float,
        default=1,
        help="minimum raw deduplicated UMI count per eligible cell",
    )
    parser.add_argument("--folds", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=128)
    parser.add_argument("--rank", type=int, default=64)
    parser.add_argument("--power-iter", type=int, default=10)
    parser.add_argument("--scale-power-iter", type=int, default=10)
    parser.add_argument("--batch-cells", type=int, default=16_384)
    parser.add_argument(
        "--data-backend",
        choices=["auto", "cuda", "pinned"],
        default="auto",
    )
    parser.add_argument("--polish-max-iter", type=int, default=32)
    parser.add_argument("--exact-inner-steps", type=int, default=32)
    parser.add_argument("--min-eq", type=float, default=5)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--nonnegative-penalty", type=float, default=1.0)
    parser.add_argument("--max-grid-expansions", type=int, default=4)
    parser.add_argument("--grid-expansion-factor", type=float, default=4.0)
    parser.add_argument(
        "--selection-rule",
        choices=[
            "minimum", "one_standard_error", "one_se_variance_retention",
        ],
        default="minimum",
    )
    parser.add_argument("--require-converged", action="store_true")
    parser.add_argument("--require-nondegenerate", action="store_true")
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--min-profile-active-fraction", type=float, default=0.9)
    parser.add_argument("--min-profile-relative-variance", type=float, default=1e-6)
    parser.add_argument("--profile-variance-retention", type=float, default=0.9)
    args = parser.parse_args()
    if args.method == "factorized":
        if not args.rank_candidate:
            raise ValueError("factorized rank CV requires --rank-candidate")
    elif not args.multiplier:
        raise ValueError(f"{args.method} CV requires --multiplier")

    if args.primer_pairs is None:
        prepared = glm_cv.prepare_alevin_glm_data(
            args.alevin_dir,
            args.salmon_ref,
            ec_design=args.design,
            regularization_target="theta",
            min_eq=args.min_eq,
        )
    else:
        prepared = glm_cv.prepare_paired_primer_glm_data(
            args.alevin_dir,
            args.salmon_ref,
            args.primer_pairs,
            ec_design=args.design,
            regularization_target="theta",
            min_eq=args.min_eq,
            min_half_umis=args.min_half_umis,
            primer_sampling_model=args.primer_sampling_model,
        )
    eligible = glm_cv.sample_cells_by_count(
        prepared.counts,
        0,
        min_count=args.min_cell_umis,
        totals=prepared.cell_umi_totals,
    )
    selected = glm_cv.sample_cells_by_count(
        prepared.counts,
        args.cells,
        seed=args.seed,
        min_count=args.min_cell_umis,
        totals=prepared.cell_umi_totals,
    )
    if not len(selected):
        raise ValueError("no cells meet the requested UMI threshold")
    counts = prepared.counts[selected]
    fold_pairs = None
    if prepared.metadata and prepared.metadata.get("paired_primers"):
        if prepared.cv_raw_counts is None:
            raise ValueError("paired-primer CV requires raw molecule counts")
        fold_pairs = glm_cv.paired_primer_count_fold_pairs(
            prepared.cv_raw_counts[selected],
            n_folds=args.folds,
            seed=args.seed,
        )
    fit_kwargs = {
        "rank": args.rank,
        "max_iter": args.max_iter,
        "min_iter": min(50, args.max_iter),
        "patience": 10,
        "tol": args.tol,
    }
    if args.method == "admm_factorized":
        fit_kwargs.update(
            adaptive_rho=True,
            rho_update_interval=10,
            rho_balance=10.0,
            rho_scale=2.0,
        )
    elif args.method == "frank_wolfe_penalized":
        fit_kwargs.update(
            max_atoms=args.max_iter,
            power_iter=args.power_iter,
            nonnegative_penalty=args.nonnegative_penalty,
        )
    elif args.method == "factorized":
        fit_kwargs.update(
            minibatch=False,
            polish_max_iter=args.polish_max_iter,
            exact_inner_steps=args.exact_inner_steps,
        )
    if args.method == "factorized":
        rank_selection_rule = (
            args.selection_rule
            if args.selection_rule in {"minimum", "one_standard_error"}
            else "one_standard_error"
        )
        report = glm_cv.cross_validate_factorized_rank_adaptive(
            counts,
            prepared.compatibility,
            args.rank_candidate,
            n_folds=args.folds,
            seed=args.seed,
            device=args.device,
            batch_cells=args.batch_cells,
            data_backend=args.data_backend,
            fit_kwargs=fit_kwargs,
            min_profile_active_fraction=args.min_profile_active_fraction,
            min_profile_relative_variance=args.min_profile_relative_variance,
            max_grid_expansions=args.max_grid_expansions,
            max_rank=args.max_rank,
            selection_rule=rank_selection_rule,
            require_converged=args.require_converged,
            require_nondegenerate=args.require_nondegenerate,
            fold_pairs=fold_pairs,
            progress_callback=lambda row: print(
                json.dumps({"event": "rank_fold_complete", **row}),
                flush=True,
            ),
        )
    else:
        report = glm_cv.cross_validate_glm_adaptive_grid(
            counts,
            prepared.compatibility,
            args.method,
            args.multiplier,
            n_folds=args.folds,
            seed=args.seed,
            device=args.device,
            batch_cells=args.batch_cells,
            data_backend=args.data_backend,
            power_iter=args.scale_power_iter,
            fit_kwargs=fit_kwargs,
            min_profile_active_fraction=args.min_profile_active_fraction,
            min_profile_relative_variance=args.min_profile_relative_variance,
            max_grid_expansions=args.max_grid_expansions,
            grid_expansion_factor=args.grid_expansion_factor,
            selection_rule=args.selection_rule,
            require_converged=args.require_converged,
            require_nondegenerate=args.require_nondegenerate,
            profile_variance_retention=args.profile_variance_retention,
            fold_pairs=fold_pairs,
        )
    full_scale = None
    if args.method != "factorized":
        full_counts = (
            counts
            if len(selected) == len(eligible)
            else prepared.counts[eligible].tocsr()
        )
        full_scale = glm_cv.hyperparameter_scale(
            full_counts,
            prepared.compatibility,
            args.method,
            device=args.device,
            batch_cells=args.batch_cells,
            power_iter=args.scale_power_iter,
            seed=args.seed,
            data_backend=args.data_backend,
        )
    report.update(
        design=args.design,
        regularization_target="theta",
        min_cell_umis=float(args.min_cell_umis),
        n_cv_cells=int(len(selected)),
        n_eligible_cells=int(len(eligible)),
        n_total_cells=int(prepared.counts.shape[0]),
        n_full_cells=int(len(eligible)),
        n_equivalence_classes=int(prepared.counts.shape[1]),
        n_transcripts=int(prepared.compatibility.shape[1]),
        response_nonzeros=int(prepared.counts.nnz),
        response_storage_bytes=glm_cv.sparse_storage_bytes(prepared.counts),
        compatibility_nonzeros=int(prepared.compatibility.nnz),
        compatibility_storage_bytes=glm_cv.sparse_storage_bytes(
            prepared.compatibility
        ),
        raw_count_nonzeros=(
            None
            if prepared.cv_raw_counts is None
            else int(prepared.cv_raw_counts.nnz)
        ),
        raw_count_storage_bytes=(
            None
            if prepared.cv_raw_counts is None
            else glm_cv.sparse_storage_bytes(prepared.cv_raw_counts)
        ),
        full_scale=full_scale,
        selected_full_value=(
            report.get("best_multiplier") * full_scale
            if report.get("best_multiplier") is not None else None
        ),
        selected_rank=report.get("best_rank", args.rank),
        data_backend=args.data_backend,
        batch_cells=int(args.batch_cells),
        max_iter=int(args.max_iter),
        tolerance=float(args.tol),
        exact_inner_steps=(
            int(args.exact_inner_steps) if args.method == "factorized" else None
        ),
        paired_primers=args.primer_pairs is not None,
        primer_pair_file=(
            str(args.primer_pairs) if args.primer_pairs is not None else None
        ),
        min_half_umis=(
            float(args.min_half_umis) if args.primer_pairs is not None else None
        ),
        primer_sampling_model=(
            args.primer_sampling_model if args.primer_pairs is not None else None
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
