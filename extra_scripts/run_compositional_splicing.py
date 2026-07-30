#!/usr/bin/env python3
"""Compare compositional and joint-block differential-splicing tests."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from extra_scripts.run_differential_splicing import (
    benjamini_hochberg,
    load_blocks,
)
from tealeaf.sc import differential


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimates", required=True, type=Path)
    parser.add_argument("--gls-tests", required=True, type=Path)
    parser.add_argument("--block-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-condition-replicates", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def read_reliable_estimates(path):
    grouped = defaultdict(list)
    with gzip.open(path, "rt") as source:
        for line in source:
            row = json.loads(line)
            if not (
                row["fit_converged"] and row["conditional_reliable"]
            ):
                continue
            grouped[(row["block_id"], row["cluster"])].append(row)
    return grouped


def condition_design(records, min_replicates):
    conditions = sorted({record["condition"] for record in records})
    counts = {
        condition: sum(
            record["condition"] == condition for record in records
        )
        for condition in conditions
    }
    conditions = [
        condition
        for condition in conditions
        if counts[condition] >= int(min_replicates)
    ]
    records = [
        record for record in records if record["condition"] in conditions
    ]
    if len(conditions) < 2 or len(records) <= len(conditions):
        return None
    index = {condition: value for value, condition in enumerate(conditions)}
    design = np.zeros((len(records), len(conditions)), dtype=float)
    design[:, 0] = 1.0
    for row, record in enumerate(records):
        column = index[record["condition"]]
        if column:
            design[row, column] = 1.0
    return records, conditions, design


def effective_path_counts(records):
    counts = []
    sizes = []
    for record in records:
        proportions = np.asarray(record["path_proportions"], dtype=float)
        size = differential.effective_multinomial_size(
            proportions,
            np.asarray(record["conditional_covariance"], dtype=float),
            maximum=record["gene_umis"],
        )
        counts.append(size * proportions)
        sizes.append(size)
    return np.asarray(counts), np.asarray(sizes)


def cauchy_combination(p_values):
    p_values = np.clip(np.asarray(p_values, dtype=float), 1e-15, 1 - 1e-15)
    statistic = np.mean(np.tan(np.pi * (0.5 - p_values)))
    return float(
        np.clip(0.5 - np.arctan(statistic) / np.pi, 0.0, 1.0)
    )


def block_equivalence(path):
    blocks = load_blocks(path, None)
    representatives = {}
    block_to_representative = {}
    members = defaultdict(list)
    for block in blocks:
        groups = tuple(sorted(
            tuple(sorted(
                transcript
                for transcript, path_index in zip(
                    block.transcripts, block.path_index
                )
                if path_index == represented_path
            ))
            for represented_path in sorted(set(block.path_index))
        ))
        key = (block.gene_id, groups)
        representative = representatives.setdefault(key, block.block_id)
        block_to_representative[block.block_id] = representative
        members[representative].append(block.block_id)
    return block_to_representative, members


def collapse_equivalent(table, block_to_representative, members, identity):
    table = table.copy()
    table["representative_block_id"] = table["block_id"].map(
        block_to_representative
    )
    table["representative_block_id"] = table[
        "representative_block_id"
    ].fillna(table["block_id"])
    table = table.sort_values("block_id").drop_duplicates(
        ["representative_block_id", *identity],
        keep="first",
    )
    table["block_id"] = table["representative_block_id"]
    table["equivalent_blocks"] = table["block_id"].map(
        lambda value: ",".join(members.get(value, [value]))
    )
    return table.drop(columns="representative_block_id")


def collapse_grouped_estimates(grouped, block_to_representative):
    result = {}
    for (block_id, cluster), records in sorted(grouped.items()):
        representative = block_to_representative.get(block_id, block_id)
        result.setdefault((representative, cluster), records)
    return result


def apply_fdr(table):
    table = table.copy()
    table["fdr"] = np.nan
    for method, positions in table.groupby("method").groups.items():
        eligible = positions
        if "converged" in table:
            eligible = positions[table.loc[positions, "converged"]]
        table.loc[eligible, "fdr"] = benjamini_hochberg(
            table.loc[eligible, "p_value"].to_numpy()
        )
    return table


def fit_compositional_tests(args, grouped):
    rows = []
    null_p_values = defaultdict(list)
    rng = np.random.default_rng(args.seed)
    for (block_id, cluster), unfiltered in grouped.items():
        prepared = condition_design(
            unfiltered, args.min_condition_replicates
        )
        if prepared is None:
            continue
        records, conditions, design = prepared
        path_counts, effective_sizes = effective_path_counts(records)
        null_design = np.ones((len(records), 1))
        methods = {
            "dirichlet_multinomial": (
                differential.dirichlet_multinomial_test
            ),
            "multinomial": differential.multinomial_test,
        }
        for method, test in methods.items():
            result = test(
                path_counts,
                null_design,
                design,
                max_iter=args.max_iter,
            )
            converged = bool(
                result["null_converged"]
                and result["alternative_converged"]
            )
            row = {
                "block_id": block_id,
                "cluster": cluster,
                "method": method,
                "fitted_model": result["model"],
                "conditions": ",".join(conditions),
                "n_samples": len(records),
                "n_paths": path_counts.shape[1],
                "median_effective_count": np.median(effective_sizes),
                "statistic": result["statistic"],
                "degrees_of_freedom": result["degrees_of_freedom"],
                "p_value": result["p_value"],
                "converged": converged,
                "valid_for_inference": (
                    method == "dirichlet_multinomial"
                ),
            }
            if method == "dirichlet_multinomial":
                row["null_concentration"] = result["null_concentration"]
            permutation_statistics = []
            if converged:
                for _ in range(int(args.permutations)):
                    permuted = design[rng.permutation(len(design))]
                    null = test(
                        path_counts,
                        null_design,
                        permuted,
                        max_iter=args.max_iter,
                    )
                    if (
                        null["null_converged"]
                        and null["alternative_converged"]
                    ):
                        permutation_statistics.append(null["statistic"])
                        null_p_values[method].append(null["p_value"])
            row["permutation_count"] = len(permutation_statistics)
            row["permutation_p_value"] = (
                (
                    1
                    + np.sum(
                        np.asarray(permutation_statistics)
                        >= result["statistic"]
                    )
                )
                / (len(permutation_statistics) + 1)
                if permutation_statistics
                else np.nan
            )
            rows.append(row)
    table = pd.DataFrame(rows)
    table["fdr"] = np.nan
    for method, positions in table.groupby("method").groups.items():
        eligible = positions[table.loc[positions, "converged"]]
        table.loc[eligible, "fdr"] = benjamini_hochberg(
            table.loc[eligible, "p_value"].to_numpy()
        )
    return table, null_p_values


def joint_block_tests(
    compositional,
    gls_path,
    block_to_representative,
    members,
):
    rows = []
    gls = pd.read_csv(gls_path, sep="\t")
    gls = gls[gls["covariance_method"] == "conditional"].copy()
    gls["method"] = "gls_acat"
    gls = gls.rename(columns={"p_value": "source_p_value"})
    gls = collapse_equivalent(
        gls,
        block_to_representative,
        members,
        identity=["cluster", "method"],
    )
    sources = [gls[[
        "block_id", "cluster", "method", "source_p_value"
    ]]]
    fitted = compositional[compositional["converged"]].copy()
    fitted["method"] = fitted["method"] + "_acat"
    fitted = fitted.rename(columns={"p_value": "source_p_value"})
    sources.append(fitted[[
        "block_id", "cluster", "method", "source_p_value"
    ]])
    source = pd.concat(sources, ignore_index=True)
    for (block_id, method), group in source.groupby(
        ["block_id", "method"], sort=False
    ):
        rows.append({
            "block_id": block_id,
            "method": method,
            "cell_type_tests": len(group),
            "p_value": cauchy_combination(group["source_p_value"]),
            "valid_for_inference": method != "multinomial_acat",
        })
    table = pd.DataFrame(rows)
    table["equivalent_blocks"] = table["block_id"].map(
        lambda value: ",".join(members.get(value, [value]))
    )
    return apply_fdr(table).sort_values(["method", "p_value"])


def joint_shared_design(records, min_replicates):
    subjects = {
        (record["condition"], record["mouse"]): record["condition"]
        for record in records
    }
    condition_counts = defaultdict(int)
    for condition in subjects.values():
        condition_counts[condition] += 1
    conditions = sorted(
        condition
        for condition, count in condition_counts.items()
        if count >= int(min_replicates)
    )
    records = [
        record for record in records if record["condition"] in conditions
    ]
    cell_types = sorted({record["cluster"] for record in records})
    if len(conditions) < 2 or not records:
        return None
    cell_type_index = {
        cell_type: index for index, cell_type in enumerate(cell_types)
    }
    null_design = np.zeros((len(records), len(cell_types)), dtype=float)
    for row, record in enumerate(records):
        null_design[row, cell_type_index[record["cluster"]]] = 1.0
    subject_keys = sorted({
        (record["condition"], record["mouse"]) for record in records
    })
    subject_index = {subject: index for index, subject in enumerate(subject_keys)}
    row_subject = np.asarray([
        subject_index[(record["condition"], record["mouse"])]
        for record in records
    ])
    subject_conditions = np.asarray([
        subject[0] for subject in subject_keys
    ])
    return (
        records,
        conditions,
        null_design,
        row_subject,
        subject_conditions,
    )


def add_condition_columns(null_design, row_conditions, conditions):
    condition_index = {
        condition: index for index, condition in enumerate(conditions)
    }
    columns = np.zeros((len(row_conditions), len(conditions) - 1))
    for row, condition in enumerate(row_conditions):
        index = condition_index.get(condition)
        if index is not None and index:
            columns[row, index - 1] = 1.0
    design = np.c_[null_design, columns]
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return None
    return design


def fit_joint_shared_tests(args, grouped):
    by_block = defaultdict(list)
    for (block_id, _), records in grouped.items():
        by_block[block_id].extend(records)
    rows = []
    null_p_values = defaultdict(list)
    rng = np.random.default_rng(args.seed + 1)
    for block_id, unfiltered in by_block.items():
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
        alternative_design = add_condition_columns(
            null_design,
            subject_conditions[row_subject],
            conditions,
        )
        if alternative_design is None:
            continue
        path_counts, effective_sizes = effective_path_counts(records)
        values = np.asarray([
            record["path_logratios"] for record in records
        ])
        covariances = np.asarray([
            record["conditional_covariance"] for record in records
        ])
        fitters = {
            "joint_shared_clustered_gls": lambda design: (
                differential.clustered_multivariate_gls_test(
                    values,
                    covariances,
                    design,
                    tested_columns=range(
                        null_design.shape[1], design.shape[1]
                    ),
                    clusters=row_subject,
                )
            ),
            "joint_shared_gls": lambda design: (
                differential.multivariate_gls_test(
                    values,
                    covariances,
                    design,
                    tested_columns=range(
                        null_design.shape[1], design.shape[1]
                    ),
                )
            ),
            "joint_shared_dirichlet_multinomial": lambda design: (
                differential.dirichlet_multinomial_test(
                    path_counts,
                    null_design,
                    design,
                    max_iter=args.max_iter,
                )
            ),
            "joint_shared_multinomial": lambda design: (
                differential.multinomial_test(
                    path_counts,
                    null_design,
                    design,
                    max_iter=args.max_iter,
                )
            ),
        }
        for method, fitter in fitters.items():
            result = fitter(alternative_design)
            converged = bool(
                result.get("null_converged", True)
                and result.get("alternative_converged", True)
            )
            permutation_statistics = []
            if converged:
                for _ in range(int(args.permutations)):
                    permuted_conditions = rng.permutation(
                        subject_conditions
                    )
                    permuted_design = add_condition_columns(
                        null_design,
                        permuted_conditions[row_subject],
                        conditions,
                    )
                    if permuted_design is None:
                        continue
                    if method == "joint_shared_clustered_gls":
                        null = differential.clustered_multivariate_gls_test(
                            values,
                            covariances,
                            permuted_design,
                            tested_columns=range(
                                null_design.shape[1],
                                permuted_design.shape[1],
                            ),
                            clusters=row_subject,
                            variance_components=(
                                result["cluster_variance"],
                                result["residual_variance"],
                            ),
                        )
                    else:
                        null = fitter(permuted_design)
                    null_converged = bool(
                        null.get("null_converged", True)
                        and null.get("alternative_converged", True)
                    )
                    if null_converged:
                        permutation_statistics.append(null["statistic"])
                        null_p_values[method].append(null["p_value"])
            row = {
                "block_id": block_id,
                "method": method,
                "conditions": ",".join(conditions),
                "cell_types": null_design.shape[1],
                "subjects": len(subject_conditions),
                "samples": len(records),
                "n_paths": path_counts.shape[1],
                "median_effective_count": np.median(effective_sizes),
                "statistic": result["statistic"],
                "degrees_of_freedom": result["degrees_of_freedom"],
                "p_value": result["p_value"],
                "converged": converged,
                "valid_for_inference": (
                    method == "joint_shared_clustered_gls"
                ),
                "permutation_count": len(permutation_statistics),
                "permutation_p_value": (
                    (
                        1
                        + np.sum(
                            np.asarray(permutation_statistics)
                            >= result["statistic"]
                        )
                    )
                    / (len(permutation_statistics) + 1)
                    if permutation_statistics
                    else np.nan
                ),
            }
            if "cluster_variance" in result:
                row["cluster_variance"] = result["cluster_variance"]
                row["residual_variance"] = result["residual_variance"]
            if "null_concentration" in result:
                row["null_concentration"] = result["null_concentration"]
            rows.append(row)
    return pd.DataFrame(rows), null_p_values


def joint_gene_tests(*tables):
    sources = []
    for table in tables:
        source = (
            table[table["converged"]].copy()
            if "converged" in table
            else table.copy()
        )
        source["gene_id"] = source["block_id"].str.rsplit(
            ":B", n=1
        ).str[0]
        sources.append(source[[
            "gene_id", "method", "p_value", "valid_for_inference"
        ]])
    source = pd.concat(sources, ignore_index=True)
    rows = []
    for (gene_id, method), group in source.groupby(
        ["gene_id", "method"], sort=False
    ):
        rows.append({
            "gene_id": gene_id,
            "method": method + "_gene_acat",
            "partition_tests": len(group),
            "p_value": cauchy_combination(group["p_value"]),
            "valid_for_inference": bool(
                group["valid_for_inference"].all()
            ),
        })
    return apply_fdr(pd.DataFrame(rows)).sort_values(
        ["method", "p_value"]
    )


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    block_to_representative, members = block_equivalence(args.block_cache)
    grouped = read_reliable_estimates(args.estimates)
    grouped = collapse_grouped_estimates(grouped, block_to_representative)
    compositional, null_p_values = fit_compositional_tests(args, grouped)
    compositional = collapse_equivalent(
        compositional,
        block_to_representative,
        members,
        identity=["cluster", "method"],
    )
    compositional = apply_fdr(compositional)
    compositional.to_csv(
        args.output_dir / "differential_compositional.tsv",
        sep="\t",
        index=False,
    )
    joint = joint_block_tests(
        compositional,
        args.gls_tests,
        block_to_representative,
        members,
    )
    joint.to_csv(
        args.output_dir / "differential_joint_block.tsv",
        sep="\t",
        index=False,
    )
    joint_shared, joint_null_p_values = fit_joint_shared_tests(args, grouped)
    joint_shared = collapse_equivalent(
        joint_shared,
        block_to_representative,
        members,
        identity=["method"],
    )
    joint_shared = apply_fdr(joint_shared).sort_values(
        ["method", "p_value"]
    )
    joint_shared.to_csv(
        args.output_dir / "differential_joint_shared.tsv",
        sep="\t",
        index=False,
    )
    joint_gene = joint_gene_tests(joint, joint_shared)
    joint_gene.to_csv(
        args.output_dir / "differential_joint_gene.tsv",
        sep="\t",
        index=False,
    )
    summary = {
        "seconds": time.perf_counter() - started,
        "block_celltype_tests": {},
        "joint_block_tests": {},
        "joint_shared_tests": {},
        "joint_gene_tests": {},
    }
    for method, group in compositional.groupby("method"):
        null = np.asarray(null_p_values[method], dtype=float)
        summary["block_celltype_tests"][method] = {
            "tests": int(len(group)),
            "converged": int(group["converged"].sum()),
            "nominal_p_lt_0.05": int(
                ((group["p_value"] < 0.05) & group["converged"]).sum()
            ),
            "fdr_0.05": int((group["fdr"] <= 0.05).sum()),
            "permutation_tests": int(len(null)),
            "permutation_p_lt_0.05": (
                float(np.mean(null < 0.05)) if len(null) else None
            ),
        }
    for method, group in joint.groupby("method"):
        summary["joint_block_tests"][method] = {
            "tests": int(len(group)),
            "nominal_p_lt_0.05": int((group["p_value"] < 0.05).sum()),
            "fdr_0.05": int((group["fdr"] <= 0.05).sum()),
        }
    for method, group in joint_shared.groupby("method"):
        null = np.asarray(joint_null_p_values[method], dtype=float)
        summary["joint_shared_tests"][method] = {
            "tests": int(len(group)),
            "converged": int(group["converged"].sum()),
            "nominal_p_lt_0.05": int(
                ((group["p_value"] < 0.05) & group["converged"]).sum()
            ),
            "fdr_0.05": int((group["fdr"] <= 0.05).sum()),
            "permutation_tests": int(len(null)),
            "permutation_p_lt_0.05": (
                float(np.mean(null < 0.05)) if len(null) else None
            ),
        }
    for method, group in joint_gene.groupby("method"):
        summary["joint_gene_tests"][method] = {
            "tests": int(len(group)),
            "nominal_p_lt_0.05": int((group["p_value"] < 0.05).sum()),
            "fdr_0.05": int((group["fdr"] <= 0.05).sum()),
        }
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
