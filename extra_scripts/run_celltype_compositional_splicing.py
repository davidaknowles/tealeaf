#!/usr/bin/env python3
"""Test paired cell-type splicing with a compositional likelihood."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
from itertools import combinations
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from extra_scripts.run_compositional_splicing import (
    apply_fdr,
    block_equivalence,
    effective_path_counts,
)
from tealeaf.sc import differential
from tealeaf.sc.ds_benchmark import filter_grouped_subject_records


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimates", required=True, type=Path)
    parser.add_argument("--block-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--subject-folds", type=Path)
    parser.add_argument("--subject-fold", type=int)
    parser.add_argument("--min-celltype-mice", type=int, default=3)
    parser.add_argument("--min-path-proportion", type=float, default=0.0)
    parser.add_argument(
        "--max-paths",
        type=int,
        help=(
            "Retain this many molecule-weighted dominant paths before "
            "applying the path-proportion filter."
        ),
    )
    parser.add_argument("--max-logratio-variance", type=float, default=np.inf)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument(
        "--celltype-test",
        choices=("omnibus", "pairwise"),
        default="omnibus",
    )
    parser.add_argument(
        "--likelihood",
        choices=("dirichlet-multinomial", "multinomial"),
        default="dirichlet-multinomial",
    )
    parser.add_argument(
        "--alternative-concentration",
        choices=("fixed", "free"),
        default="fixed",
        help=(
            "Fix alternative kappa at its null estimate or re-estimate it "
            "under the alternative."
        ),
    )
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def read_identifiable_estimates(
    path,
    block_to_representative,
    min_path_proportion,
    max_logratio_variance,
    max_paths=None,
):
    grouped = defaultdict(list)
    partition_source = {}
    with gzip.open(path, "rt") as source:
        for line in source:
            row = json.loads(line)
            if not (
                row["fit_converged"] and row["conditional_identifiable"]
            ):
                continue
            representative = block_to_representative.get(
                row["block_id"], row["block_id"]
            )
            source_block = partition_source.setdefault(
                representative, row["block_id"]
            )
            if row["block_id"] != source_block:
                continue
            grouped[representative].append(row)
    projected = {}
    for block_id, records in grouped.items():
        n_paths = len(records[0]["path_proportions"])
        selected = np.arange(n_paths)
        if max_paths is not None and n_paths > int(max_paths):
            abundance = sum(
                float(record["gene_umis"])
                * np.asarray(record["path_proportions"], dtype=float)
                for record in records
            )
            selected = np.sort(
                np.argsort(abundance)[-int(max_paths) :]
            )
        retained = []
        for record in records:
            proportions = np.asarray(
                record["path_proportions"], dtype=float
            )
            covariance = np.asarray(
                record["conditional_covariance"], dtype=float
            )
            if len(selected) != n_paths:
                proportions, covariance = (
                    differential.project_composition(
                        proportions, covariance, selected
                    )
                )
            if (
                proportions.min() < float(min_path_proportion)
                or np.linalg.eigvalsh(covariance).max()
                > float(max_logratio_variance)
            ):
                continue
            record = dict(record)
            record["path_proportions"] = proportions.tolist()
            record["conditional_covariance"] = covariance.tolist()
            record["selected_path_indices"] = selected.tolist()
            retained.append(record)
        if retained:
            projected[block_id] = retained
    return projected


def paired_celltype_design(records, min_mice):
    counts = Counter(
        (record["cluster"], record["condition"], record["mouse"])
        for record in records
    )
    cell_type_subjects = Counter(
        cell_type for cell_type, _, _ in counts
    )
    cell_types = sorted(
        cell_type
        for cell_type, count in cell_type_subjects.items()
        if count >= int(min_mice)
    )
    records = [
        record for record in records if record["cluster"] in cell_types
    ]
    if len(cell_types) < 2:
        return None
    subjects = sorted({
        (record["condition"], record["mouse"]) for record in records
    })
    subject_index = {
        subject: index for index, subject in enumerate(subjects)
    }
    cell_type_index = {
        cell_type: index for index, cell_type in enumerate(cell_types)
    }
    null_design = np.zeros((len(records), len(subjects)), dtype=float)
    null_design[:, 0] = 1.0
    alternative_design = np.zeros(
        (len(records), len(subjects) + len(cell_types) - 1),
        dtype=float,
    )
    alternative_design[:, : len(subjects)] = null_design
    row_subject = np.empty(len(records), dtype=np.int64)
    for row, record in enumerate(records):
        subject = (record["condition"], record["mouse"])
        subject_column = subject_index[subject]
        row_subject[row] = subject_column
        if subject_column:
            null_design[row, subject_column] = 1.0
            alternative_design[row, subject_column] = 1.0
        cell_type_column = cell_type_index[record["cluster"]]
        if cell_type_column:
            alternative_design[
                row, len(subjects) + cell_type_column - 1
            ] = 1.0
    if (
        len(records) <= alternative_design.shape[1]
        or np.linalg.matrix_rank(alternative_design)
        < alternative_design.shape[1]
    ):
        return None
    return (
        records,
        cell_types,
        subjects,
        row_subject,
        null_design,
        alternative_design,
    )


def pairwise_celltype_designs(records, min_mice):
    by_celltype_subject = {
        (record["cluster"], record["condition"], record["mouse"]): record
        for record in records
    }
    cell_types = sorted({record["cluster"] for record in records})
    designs = []
    for first, second in combinations(cell_types, 2):
        first_subjects = {
            subject[1:]
            for subject in by_celltype_subject
            if subject[0] == first
        }
        second_subjects = {
            subject[1:]
            for subject in by_celltype_subject
            if subject[0] == second
        }
        subjects = sorted(first_subjects & second_subjects)
        if len(subjects) < int(min_mice):
            continue
        paired_records = [
            by_celltype_subject[(cell_type, *subject)]
            for subject in subjects
            for cell_type in (first, second)
        ]
        prepared = paired_celltype_design(paired_records, min_mice)
        if prepared is not None:
            designs.append((f"{first}_vs_{second}", prepared))
    return designs


def fit_tests(args, grouped, members):
    rows = []
    null_records = []
    rng = np.random.default_rng(args.seed)
    items = sorted(grouped.items())
    if args.max_blocks is not None:
        items = items[: int(args.max_blocks)]
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("shard index must be between zero and shard count")
    items = items[args.shard_index :: args.shard_count]
    for block_id, unfiltered in items:
        if args.celltype_test == "pairwise":
            designs = pairwise_celltype_designs(
                unfiltered, args.min_celltype_mice
            )
        else:
            prepared = paired_celltype_design(
                unfiltered, args.min_celltype_mice
            )
            designs = [("omnibus", prepared)] if prepared is not None else []
        for contrast, prepared in designs:
            (
                records,
                cell_types,
                subjects,
                row_subject,
                null_design,
                alternative_design,
            ) = prepared
            path_counts, effective_sizes = effective_path_counts(records)
            test = (
                differential.dirichlet_multinomial_test
                if args.likelihood == "dirichlet-multinomial"
                else differential.multinomial_test
            )
            test_kwargs = {}
            if args.likelihood == "dirichlet-multinomial":
                test_kwargs["fix_null_concentration"] = (
                    args.alternative_concentration == "fixed"
                )
            result = test(
                path_counts,
                null_design,
                alternative_design,
                max_iter=args.max_iter,
                **test_kwargs,
            )
            converged = bool(
                result["null_converged"]
                and result["alternative_converged"]
            )
            cell_type_effects = result["alternative_coefficients"][
                null_design.shape[1] :
            ]
            effect_norms = np.linalg.norm(cell_type_effects, axis=1)
            permutation_statistics = []
            test_id = f"{block_id}|{contrast}"
            if converged:
                subject_rows = [
                    np.flatnonzero(row_subject == subject)
                    for subject in range(len(subjects))
                ]
                tested = np.arange(
                    null_design.shape[1], alternative_design.shape[1]
                )
                for _ in range(int(args.permutations)):
                    permuted_design = alternative_design.copy()
                    for positions in subject_rows:
                        shuffled = rng.permutation(positions)
                        permuted_design[
                            positions[:, None], tested
                        ] = alternative_design[shuffled[:, None], tested]
                    null = test(
                        path_counts,
                        null_design,
                        permuted_design,
                        max_iter=args.max_iter,
                        fitted_null=result,
                        **test_kwargs,
                    )
                    if (
                        null["null_converged"]
                        and null["alternative_converged"]
                    ):
                        permutation_statistics.append(null["statistic"])
                        null_records.append({
                            "test_id": test_id,
                            "block_id": block_id,
                            "degrees_of_freedom": result[
                                "degrees_of_freedom"
                            ],
                            "p_value": null["p_value"],
                            "statistic": null["statistic"],
                        })
            permutation_p_value = (
                differential.permutation_rank_p_value(
                    result["statistic"], permutation_statistics
                )
                if permutation_statistics
                else np.nan
            )
            rows.append({
                "test_id": test_id,
                "block_id": block_id,
                "equivalent_blocks": ",".join(
                    members.get(block_id, [block_id])
                ),
                "contrast": contrast,
                "method": f"paired_{args.likelihood.replace('-', '_')}",
                "fitted_model": result["model"],
                "cell_types": ",".join(cell_types),
                "n_cell_types": len(cell_types),
                "n_mice": len(subjects),
                "n_samples": len(records),
                "n_paths": path_counts.shape[1],
                "selected_path_indices": ",".join(
                    map(str, records[0]["selected_path_indices"])
                ),
                "median_effective_count": np.median(effective_sizes),
                "statistic": result["statistic"],
                "degrees_of_freedom": result["degrees_of_freedom"],
                "asymptotic_p_value": result["p_value"],
                "p_value": (
                    permutation_p_value
                    if permutation_statistics
                    else result["p_value"]
                ),
                "converged": converged,
                "null_concentration": result.get(
                    "null_concentration", np.nan
                ),
                "alternative_concentration": result.get(
                    "alternative_concentration", np.nan
                ),
                "celltype_logit_effects": json.dumps(
                    cell_type_effects.tolist()
                ),
                "max_celltype_logit_effect": (
                    float(effect_norms.max())
                    if len(effect_norms)
                    else 0.0
                ),
                "permutation_count": len(permutation_statistics),
                "permutation_p_value": permutation_p_value,
                "valid_for_inference": True,
            })
    table = apply_fdr(pd.DataFrame(rows))
    return table.sort_values("p_value"), pd.DataFrame(null_records)


def main():
    args = parse_args()
    if args.max_paths is not None and args.max_paths < 2:
        raise ValueError("--max-paths must be at least two")
    if (
        args.likelihood != "dirichlet-multinomial"
        and args.alternative_concentration != "fixed"
    ):
        raise ValueError(
            "--alternative-concentration applies only to "
            "Dirichlet-multinomial fits"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    block_to_representative, members = block_equivalence(args.block_cache)
    grouped = read_identifiable_estimates(
        args.estimates,
        block_to_representative,
        args.min_path_proportion,
        args.max_logratio_variance,
        args.max_paths,
    )
    selected_subjects = None
    if args.subject_folds is not None:
        if args.subject_fold is None:
            raise ValueError("--subject-fold is required with --subject-folds")
        folds = pd.read_csv(
            args.subject_folds, sep="\t", dtype={"subject": str}
        )
        selected_subjects = set(
            folds.loc[folds.fold.eq(args.subject_fold), "subject"]
        )
        grouped = filter_grouped_subject_records(grouped, selected_subjects)
    elif args.subject_fold is not None:
        raise ValueError("--subject-folds is required with --subject-fold")
    table, null_table = fit_tests(args, grouped, members)
    table.to_csv(
        args.output_dir / "differential_cell_type_compositional.tsv",
        sep="\t",
        index=False,
    )
    null_table.to_csv(
        args.output_dir / "permutation_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summary = {
        "seconds": time.perf_counter() - started,
        "min_path_proportion": args.min_path_proportion,
        "max_logratio_variance": args.max_logratio_variance,
        "max_paths": args.max_paths,
        "celltype_test": args.celltype_test,
        "likelihood": args.likelihood,
        "alternative_concentration": args.alternative_concentration,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "candidate_partitions": len(grouped),
        "subject_fold": args.subject_fold,
        "subjects": (
            len(selected_subjects) if selected_subjects is not None else None
        ),
        "tests": int(len(table)),
        "converged": int(table["converged"].sum()),
        "nominal_p_lt_0.05": int(
            ((table["p_value"] < 0.05) & table["converged"]).sum()
        ),
        "fdr_0.05": int((table["fdr"] <= 0.05).sum()),
        "permutation_tests": int(len(null_table)),
        "permutation_p_lt_0.05": (
            float(np.mean(null_table["p_value"] < 0.05))
            if len(null_table)
            else None
        ),
    }
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
