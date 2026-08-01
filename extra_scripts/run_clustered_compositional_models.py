#!/usr/bin/env python3
"""Compare clustered multinomial GLMM, GEE, and existing splicing tests."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy.special

from extra_scripts.run_compositional_splicing import (
    add_condition_columns,
    apply_fdr,
    block_equivalence,
    collapse_grouped_estimates,
    effective_path_counts,
    joint_shared_design,
    read_reliable_estimates,
)
from tealeaf.sc import differential


METHODS = (
    "gee_exchangeable",
    "glmm_laplace",
    "clustered_gls",
    "dirichlet_multinomial",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimates", required=True, type=Path)
    parser.add_argument("--block-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-condition-replicates", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--power-replicates", type=int, default=20)
    parser.add_argument(
        "--effect-sizes", nargs="+", type=float, default=(0.25, 0.5)
    )
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def fit_method(
    method,
    path_counts,
    values,
    covariances,
    null_design,
    alternative_design,
    row_subject,
    *,
    max_iter,
    fitted_null=None,
):
    tested = range(null_design.shape[1], alternative_design.shape[1])
    if method == "gee_exchangeable":
        result = differential.multinomial_gee_test(
            path_counts,
            alternative_design,
            tested,
            row_subject,
            max_iter=max_iter,
            initial=(
                fitted_null["coefficients"]
                if fitted_null is not None
                else None
            ),
        )
        result["fit_converged"] = result["converged"]
    elif method == "glmm_laplace":
        integer_counts = differential.integerize_compositional_counts(
            path_counts
        )
        result = differential.multinomial_glmm_test(
            integer_counts,
            null_design,
            alternative_design,
            row_subject,
            max_iter=max_iter,
            fitted_null=fitted_null,
        )
        result["fit_converged"] = bool(
            result["null_converged"]
            and result["alternative_converged"]
        )
    elif method == "clustered_gls":
        result = differential.clustered_multivariate_gls_test(
            values,
            covariances,
            alternative_design,
            tested,
            row_subject,
            variance_components=(
                (
                    fitted_null["cluster_variance"],
                    fitted_null["residual_variance"],
                )
                if fitted_null is not None
                else None
            ),
        )
        result["fit_converged"] = True
    elif method == "dirichlet_multinomial":
        result = differential.dirichlet_multinomial_test(
            path_counts,
            null_design,
            alternative_design,
            max_iter=max_iter,
            fitted_null=fitted_null,
        )
        result["fit_converged"] = bool(
            result["null_converged"]
            and result["alternative_converged"]
        )
    else:  # pragma: no cover - guarded by argparse
        raise ValueError(f"unknown method {method}")
    return result


def shifted_composition(values, effect, condition, direction):
    shifted = values + float(effect) * condition[:, None] * direction[None, :]
    basis = differential.helmert_basis(values.shape[1] + 1)
    centered_logits = shifted @ basis.T
    return scipy.special.softmax(centered_logits, axis=1), shifted


def result_row(block_id, method, result, metadata):
    row = {
        "block_id": block_id,
        "method": method,
        **metadata,
        "statistic": result.get("statistic", np.nan),
        "degrees_of_freedom": result.get("degrees_of_freedom", 0),
        "p_value": result.get("p_value", np.nan),
        "converged": bool(result.get("fit_converged", False)),
    }
    for name in (
        "working_correlation",
        "scale",
        "null_random_effect_sd",
        "alternative_random_effect_sd",
        "cluster_variance",
        "residual_variance",
        "null_concentration",
        "alternative_concentration",
        "null_gradient_norm",
        "alternative_gradient_norm",
        "null_mode_score_norm",
        "alternative_mode_score_norm",
    ):
        if name in result:
            row[name] = result[name]
    return row


def evaluate(args, grouped):
    by_block = defaultdict(list)
    for (block_id, _), records in grouped.items():
        by_block[block_id].extend(records)
    items = sorted(by_block.items())
    if args.max_blocks is not None:
        items = items[: int(args.max_blocks)]
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("shard index must be between zero and shard count")
    items = items[args.shard_index :: args.shard_count]
    rng = np.random.default_rng(args.seed + args.shard_index)
    observed_rows = []
    null_rows = []
    power_rows = []
    failures = []
    eligible_blocks = 0
    for block_number, (block_id, unfiltered) in enumerate(items):
        prepared = joint_shared_design(
            unfiltered, args.min_condition_replicates
        )
        if prepared is None:
            continue
        (
            records,
            conditions,
            null_design,
            row_subject,
            subject_conditions,
        ) = prepared
        row_conditions = subject_conditions[row_subject]
        alternative_design = add_condition_columns(
            null_design, row_conditions, conditions
        )
        if alternative_design is None:
            continue
        eligible_blocks += 1
        path_counts, effective_sizes = effective_path_counts(records)
        values = np.asarray([
            record["path_logratios"] for record in records
        ], dtype=float)
        covariances = np.asarray([
            record["conditional_covariance"] for record in records
        ], dtype=float)
        metadata = {
            "conditions": ",".join(conditions),
            "n_conditions": len(conditions),
            "n_cell_types": null_design.shape[1],
            "n_subjects": len(subject_conditions),
            "n_samples": len(records),
            "n_paths": path_counts.shape[1],
            "median_effective_count": float(np.median(effective_sizes)),
        }
        observed_results = {}
        for method in args.methods:
            try:
                result = fit_method(
                    method,
                    path_counts,
                    values,
                    covariances,
                    null_design,
                    alternative_design,
                    row_subject,
                    max_iter=args.max_iter,
                )
                observed_rows.append(
                    result_row(block_id, method, result, metadata)
                )
                observed_results[method] = result
            except Exception as exc:
                failures.append({
                    "block_id": block_id,
                    "method": method,
                    "stage": "observed",
                    "error": repr(exc),
                })
        permutations = [
            rng.permutation(subject_conditions)
            for _ in range(int(args.permutations))
        ]
        for replicate, permuted_subject_conditions in enumerate(permutations):
            permuted_design = add_condition_columns(
                null_design,
                permuted_subject_conditions[row_subject],
                conditions,
            )
            if permuted_design is None:
                continue
            for method in args.methods:
                try:
                    result = fit_method(
                        method,
                        path_counts,
                        values,
                        covariances,
                        null_design,
                        permuted_design,
                        row_subject,
                        max_iter=args.max_iter,
                        fitted_null=observed_results.get(method),
                    )
                    null_rows.append(result_row(
                        block_id,
                        method,
                        result,
                        {**metadata, "replicate": replicate},
                    ))
                except Exception as exc:
                    failures.append({
                        "block_id": block_id,
                        "method": method,
                        "stage": f"null_{replicate}",
                        "error": repr(exc),
                    })
        if len(conditions) != 2 or not args.power_replicates:
            if eligible_blocks % 5 == 0:
                print(json.dumps({
                    "blocks_scanned": block_number + 1,
                    "eligible_blocks": eligible_blocks,
                    "observed_fits": len(observed_rows),
                    "null_fits": len(null_rows),
                    "power_fits": len(power_rows),
                    "failures": len(failures),
                }), flush=True)
            continue
        direction = rng.normal(size=values.shape[1])
        direction /= np.linalg.norm(direction)
        for replicate in range(int(args.power_replicates)):
            permuted_subject_conditions = rng.permutation(subject_conditions)
            power_design = add_condition_columns(
                null_design,
                permuted_subject_conditions[row_subject],
                conditions,
            )
            if power_design is None:
                continue
            condition = power_design[:, -1]
            for effect in args.effect_sizes:
                shifted_proportions, shifted_values = shifted_composition(
                    values, effect, condition, direction
                )
                shifted_counts = effective_sizes[:, None] * shifted_proportions
                for method in args.methods:
                    try:
                        result = fit_method(
                            method,
                            shifted_counts,
                            shifted_values,
                            covariances,
                            null_design,
                            power_design,
                            row_subject,
                            max_iter=args.max_iter,
                            fitted_null=(
                                observed_results.get(method)
                                if method == "gee_exchangeable"
                                else None
                            ),
                        )
                        power_rows.append(result_row(
                            block_id,
                            method,
                            result,
                            {
                                **metadata,
                                "replicate": replicate,
                                "effect_size": effect,
                            },
                        ))
                    except Exception as exc:
                        failures.append({
                            "block_id": block_id,
                            "method": method,
                            "stage": f"power_{effect}_{replicate}",
                            "error": repr(exc),
                        })
        if eligible_blocks % 5 == 0:
            print(json.dumps({
                "blocks_scanned": block_number + 1,
                "eligible_blocks": eligible_blocks,
                "observed_fits": len(observed_rows),
                "null_fits": len(null_rows),
                "power_fits": len(power_rows),
                "failures": len(failures),
            }), flush=True)
    return (
        pd.DataFrame(observed_rows),
        pd.DataFrame(null_rows),
        pd.DataFrame(power_rows),
        failures,
    )


def summarize(observed, null, power):
    observed = apply_fdr(observed) if len(observed) else observed
    summary = {"observed": {}, "null": {}, "power": {}}
    null_statistics = {}
    if len(null):
        for (block_id, method), group in null.groupby(
            ["block_id", "method"]
        ):
            valid = group["converged"] & np.isfinite(group["p_value"])
            statistics = group.loc[valid, "statistic"].to_numpy()
            if len(statistics):
                null_statistics[(block_id, method)] = statistics
        for method, group in null.groupby("method"):
            valid = group["converged"] & np.isfinite(group["p_value"])
            values = group.loc[valid, "p_value"].to_numpy()
            summary["null"][method] = {
                "fits": int(len(group)),
                "converged": int(valid.sum()),
                "asymptotic_rejection_0.05": (
                    float(np.mean(values < 0.05)) if len(values) else None
                ),
                "pooled_p05_quantile": (
                    float(np.quantile(values, 0.05))
                    if len(values)
                    else None
                ),
            }

    def permutation_p_values(table):
        values = []
        for row in table.itertuples():
            statistics = null_statistics.get((row.block_id, row.method))
            if statistics is None or not np.isfinite(row.statistic):
                values.append(np.nan)
            else:
                values.append(differential.permutation_rank_p_value(
                    row.statistic, statistics
                ))
        return values

    if len(observed):
        observed = observed.copy()
        observed["permutation_p_value"] = permutation_p_values(observed)
        permutation_table = observed.copy()
        permutation_table["p_value"] = observed["permutation_p_value"]
        permutation_table["converged"] = (
            observed["converged"]
            & np.isfinite(observed["permutation_p_value"])
        )
        observed["permutation_fdr"] = apply_fdr(permutation_table)["fdr"]
        for method, group in observed.groupby("method"):
            valid = group["converged"] & np.isfinite(group["p_value"])
            permutation_valid = group["converged"] & np.isfinite(
                group["permutation_p_value"]
            )
            summary["observed"][method] = {
                "tests": int(len(group)),
                "converged": int(valid.sum()),
                "nominal_p_lt_0.05": int(
                    (group.loc[valid, "p_value"] < 0.05).sum()
                ),
                "fdr_0.05": int((group["fdr"] <= 0.05).sum()),
                "permutation_p_lt_0.05": int((
                    group.loc[
                        permutation_valid, "permutation_p_value"
                    ] <= 0.05
                ).sum()),
                "permutation_fdr_0.05": int((
                    group["permutation_fdr"] <= 0.05
                ).sum()),
            }
    if len(power):
        power = power.copy()
        power["calibrated_p_value"] = permutation_p_values(power)
        power["calibrated_reject"] = power["calibrated_p_value"] <= 0.05
        for (method, effect), group in power.groupby(
            ["method", "effect_size"]
        ):
            valid = group["converged"] & np.isfinite(group["p_value"])
            summary["power"].setdefault(method, {})[str(effect)] = {
                "fits": int(len(group)),
                "converged": int(valid.sum()),
                "asymptotic_power": (
                    float(np.mean(group.loc[valid, "p_value"] < 0.05))
                    if valid.any()
                    else None
                ),
                "null_calibrated_power": (
                    float(np.mean(group.loc[
                        valid & np.isfinite(group["calibrated_p_value"]),
                        "calibrated_reject",
                    ]))
                    if (
                        valid & np.isfinite(group["calibrated_p_value"])
                    ).any()
                    else None
                ),
            }
    return observed, power, summary


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    block_to_representative, _ = block_equivalence(args.block_cache)
    grouped = read_reliable_estimates(args.estimates)
    grouped = collapse_grouped_estimates(grouped, block_to_representative)
    observed, null, power, failures = evaluate(args, grouped)
    observed, power, summary = summarize(observed, null, power)
    observed.to_csv(
        args.output_dir / "clustered_model_observed.tsv",
        sep="\t",
        index=False,
    )
    null.to_csv(
        args.output_dir / "clustered_model_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    power.to_csv(
        args.output_dir / "clustered_model_power.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    with gzip.open(args.output_dir / "failures.jsonl.gz", "wt") as output:
        for row in failures:
            output.write(json.dumps(row) + "\n")
    summary.update({
        "seconds": time.perf_counter() - started,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "permutations": args.permutations,
        "power_replicates": args.power_replicates,
        "effect_sizes": args.effect_sizes,
        "failures": len(failures),
    })
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
