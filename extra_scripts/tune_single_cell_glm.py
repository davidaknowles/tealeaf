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
    parser.add_argument("--design", required=True, choices=["binary", "weighted"])
    parser.add_argument(
        "--method",
        required=True,
        choices=["admm_factorized", "frank_wolfe_penalized"],
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--multiplier", action="append", required=True, type=float)
    parser.add_argument("--cells", type=int, default=5_000)
    parser.add_argument("--folds", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=128)
    parser.add_argument("--rank", type=int, default=64)
    parser.add_argument("--power-iter", type=int, default=10)
    parser.add_argument("--scale-power-iter", type=int, default=10)
    parser.add_argument("--batch-cells", type=int, default=2_048)
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

    prepared = glm_cv.prepare_alevin_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        ec_design=args.design,
        regularization_target="theta",
        min_eq=args.min_eq,
    )
    selected = glm_cv.sample_nonempty_cells(
        prepared.counts, args.cells, seed=args.seed
    )
    counts = prepared.counts[selected]
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
    else:
        fit_kwargs.update(
            max_atoms=args.max_iter,
            power_iter=args.power_iter,
            nonnegative_penalty=args.nonnegative_penalty,
        )
    report = glm_cv.cross_validate_glm_adaptive_grid(
        counts,
        prepared.compatibility,
        args.method,
        args.multiplier,
        n_folds=args.folds,
        seed=args.seed,
        device=args.device,
        batch_cells=args.batch_cells,
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
    )
    full_scale = glm_cv.hyperparameter_scale(
        prepared.counts,
        prepared.compatibility,
        args.method,
        device=args.device,
        batch_cells=args.batch_cells,
        power_iter=args.scale_power_iter,
        seed=args.seed,
    )
    report.update(
        design=args.design,
        regularization_target="theta",
        n_cv_cells=int(len(selected)),
        n_full_cells=int(prepared.counts.shape[0]),
        n_equivalence_classes=int(prepared.counts.shape[1]),
        n_transcripts=int(prepared.compatibility.shape[1]),
        full_scale=full_scale,
        selected_full_value=(
            report["best_multiplier"] * full_scale
            if report["best_multiplier"] is not None else None
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
