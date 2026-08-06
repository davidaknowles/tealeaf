#!/usr/bin/env python3
"""Plan and run external junction-based differential-splicing comparators."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from tealeaf.sc.junction_benchmark import (
    JunctionBundle,
    normalize_pvalue_table,
    plan_pairwise_contrasts,
    simes_omnibus,
    stratified_subject_folds,
)


def load_contrast(path: Path, index: int) -> dict:
    contrasts = json.loads(path.read_text())
    if index < 0 or index >= len(contrasts):
        raise IndexError(f"contrast index {index} outside 0..{len(contrasts)-1}")
    return contrasts[index]


def command_plan(args):
    bundle = JunctionBundle.load(args.bundle)
    samples = bundle.samples
    if args.subject_folds is not None:
        if args.subject_fold is None:
            raise ValueError("--subject-fold is required with --subject-folds")
        folds = pd.read_csv(args.subject_folds, sep="\t", dtype={"subject": str})
        selected = set(folds.loc[folds.fold.eq(args.subject_fold), "subject"])
        samples = samples.loc[samples.subject.astype(str).isin(selected)]
    contrasts = plan_pairwise_contrasts(samples, args.min_subjects)
    if args.effect:
        contrasts = [row for row in contrasts if row["effect"] in args.effect]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(contrasts, indent=2) + "\n")
    summary = pd.DataFrame(contrasts).groupby("effect").size().to_dict() if contrasts else {}
    print(f"contrasts={len(contrasts)} by_effect={summary}")


def command_split_subjects(args):
    bundle = JunctionBundle.load(args.bundle)
    folds = stratified_subject_folds(bundle.samples, args.folds, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    folds.to_csv(args.output, sep="\t", index=False)
    print(f"subjects={len(folds)} folds={args.folds}")


def scquint_adata(bundle):
    import anndata

    var = bundle.junctions.copy()
    var.index = (
        var.chromosome.astype(str)
        + ":"
        + var.start.astype(str)
        + "-"
        + var.end.astype(str)
        + ":"
        + var.strand.astype(str)
    )
    var.index = pd.Index(var.index).map(str)
    return anndata.AnnData(X=bundle.counts, obs=bundle.samples.set_index("sample_id"), var=var)


def command_scquint(args):
    if args.scquint_root is not None:
        sys.path.insert(0, str(args.scquint_root))
    import joblib
    import scquint.differential_splicing as scquint_ds
    import torch

    original_regression = scquint_ds.run_regression
    original_group_filter = scquint_ds.filter_min_cells_per_intron_group

    def conservative_regression(values):
        try:
            return original_regression(values)
        except Exception as error:
            intron_group, counts, *_ = values
            group = pd.DataFrame(
                {
                    "intron_group": [intron_group],
                    "p_value": [1.0],
                    "optimization_failed": [True],
                    "failure": [str(error)],
                }
            )
            introns = pd.DataFrame(
                {"psi_a": np.nan, "psi_b": np.nan}, index=np.arange(counts.shape[1])
            )
            return group, introns

    scquint_ds.run_regression = conservative_regression

    def allow_empty_group_filter(adata, *filter_args, **filter_kwargs):
        if adata.shape[1] == 0:
            return adata
        return original_group_filter(adata, *filter_args, **filter_kwargs)

    scquint_ds.filter_min_cells_per_intron_group = allow_empty_group_filter
    torch.set_num_threads(1)

    bundle = JunctionBundle.load(args.bundle)
    contrast = load_contrast(args.contrasts, args.contrast_index)
    adata = scquint_adata(bundle)
    if "intron_group" not in adata.var:
        raise ValueError("bundle has no scQuint annotation; aggregate it with --gtf")
    keep = adata.var.intron_group.notna() & (adata.var.intron_group_size >= 2)
    adata = adata[:, keep].copy()
    lookup = {sample: index for index, sample in enumerate(adata.obs_names)}
    first = np.asarray([lookup[sample] for sample in contrast["samples_a"] if sample in lookup])
    second = np.asarray([lookup[sample] for sample in contrast["samples_b"] if sample in lookup])
    with joblib.parallel_backend("threading", n_jobs=args.jobs):
        groups, introns = scquint_ds.run_differential_splicing(
            adata,
            first,
            second,
            n_jobs=args.jobs,
            min_cells_per_intron_group=args.min_samples,
            min_total_cells_per_intron=args.min_samples,
            min_global_proportion=args.min_global_proportion,
        )
    if groups.empty:
        output = pd.DataFrame()
    else:
        output = normalize_pvalue_table(
            groups.reset_index(),
            method="scQuint",
            contrast=contrast,
            feature_column="intron_group",
            pvalue_column="p_value",
        )
        for column in ("gene_id", "gene_name"):
            if column in groups:
                output[column] = groups.reset_index()[column].to_numpy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


def command_leafcutter_metadata(args):
    bundle = JunctionBundle.load(args.bundle)
    contrast = load_contrast(args.contrasts, args.contrast_index)
    rows = []
    subject_lookup = bundle.samples.set_index("sample_id").subject.to_dict()
    for group, samples in ((contrast["level_a"], contrast["samples_a"]), (contrast["level_b"], contrast["samples_b"])):
        for sample in samples:
            row = [sample, group]
            if contrast["effect"] == "cell_type":
                row.append("subject_" + str(subject_lookup[sample]))
            rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", header=False, index=False)


def command_normalize_leafcutter(args):
    contrast = load_contrast(args.contrasts, args.contrast_index)
    table = pd.read_csv(args.input, sep="\t")
    pvalue = next(name for name in ("p", "p_value", "p.adjust") if name in table)
    feature = next(name for name in ("cluster", "intron") if name in table)
    table = table.loc[pd.to_numeric(table[pvalue], errors="coerce").notna()].copy()
    if table.empty:
        output = pd.DataFrame()
    else:
        output = normalize_pvalue_table(
            table,
            method="LeafCutter",
            contrast=contrast,
            feature_column=feature,
            pvalue_column=pvalue,
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


def command_normalize_majiq(args):
    contrast = load_contrast(args.contrasts, args.contrast_index)
    table = pd.read_csv(args.input, sep="\t", comment="#")
    mode = getattr(args, "majiq_pvalue_mode", "raw")
    preferred = {
        "raw": (
            "mannwhitneyu-raw_pvalue",
            "raw_pvalue",
        ),
        "posterior_q95": (
            "mannwhitneyu-approximate_pvalue_quantiles_0.950",
            "approximate_pvalue_quantiles_0.950",
        ),
    }[mode]
    pvalue = next((name for name in preferred if name in table), None)
    if pvalue is None:
        raise ValueError(
            f"MAJIQ table has no {mode} p-value column: {list(table)}"
        )
    if "lsv_id" in table:
        feature_columns = [
            name for name in ("lsv_id", "ec_idx") if name in table
        ]
    else:
        feature_columns = [
            name
            for name in (
                "gene_id",
                "seqid",
                "strand",
                "event_type",
                "ref_exon_start",
                "ref_exon_end",
                "start",
                "end",
                "is_intron",
                "other_exon_start",
                "other_exon_end",
            )
            if name in table
        ]
    if not feature_columns or feature_columns == ["gene_id"]:
        table = table.reset_index(names="majiq_row")
        feature_columns = ["majiq_row"]
    if table.empty:
        table["majiq_feature"] = pd.Series(dtype=str)
    else:
        table["majiq_feature"] = table[feature_columns].astype(str).agg(
            ":".join, axis=1
        )
    duplicated = table["majiq_feature"].duplicated(keep=False)
    if duplicated.any():
        table.loc[duplicated, "majiq_feature"] += ":" + table.loc[
            duplicated
        ].groupby("majiq_feature").cumcount().astype(str)
    output = normalize_pvalue_table(
        table,
        method="MAJIQ Heterogen",
        contrast=contrast,
        feature_column="majiq_feature",
        pvalue_column=pvalue,
    )
    output["majiq_pvalue_column"] = pvalue
    for column in ("gene_id", "gene_name"):
        if column in table:
            output[column] = table[column].to_numpy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


def command_normalize_majiq_directory(args):
    """Renormalize existing MAJIQ tables without rerunning Heterogen."""
    contrasts = json.loads(args.contrasts.read_text())
    contrast_index = {
        str(contrast["contrast_id"]): index
        for index, contrast in enumerate(contrasts)
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    completed = 0
    for path in sorted(args.input_dir.glob("*.tsv")):
        index = contrast_index.get(path.stem)
        if index is None:
            raise ValueError(f"MAJIQ table {path} has no planned contrast")
        command_normalize_majiq(
            argparse.Namespace(
                contrasts=args.contrasts,
                contrast_index=index,
                input=path,
                output=args.output_dir / path.name,
                majiq_pvalue_mode=args.majiq_pvalue_mode,
            )
        )
        completed += 1
    print(f"normalized_majiq_tables={completed}")


def command_summarize(args):
    paths = list(args.input or [])
    for directory in args.input_dir or []:
        paths.extend(sorted(directory.glob("*.tsv")))
    tables = []
    for path in paths:
        if not path.stat().st_size:
            continue
        try:
            table = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if not table.empty:
            tables.append(table)
    combined = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False)
    if combined.empty:
        pd.DataFrame().to_csv(args.omnibus_output, sep="\t", index=False)
        pd.DataFrame(
            columns=["method", "effect", "tests", "discoveries"]
        ).to_csv(args.summary_output, sep="\t", index=False)
        return
    cell_type = combined[combined.effect == "cell_type"]
    omnibus = simes_omnibus(cell_type)
    omnibus.to_csv(args.omnibus_output, sep="\t", index=False)
    condition_summary = (
        combined.loc[combined.effect == "condition"]
        .groupby(["method", "effect"])
        .agg(tests=("feature_id", "size"), discoveries=("significant", "sum"))
        .reset_index()
    )
    cell_type_summary = (
        omnibus.groupby("method")
        .agg(tests=("feature_id", "size"), discoveries=("significant", "sum"))
        .reset_index()
    )
    cell_type_summary.insert(1, "effect", "cell_type")
    summary = pd.concat(
        [cell_type_summary, condition_summary], ignore_index=True
    ).sort_values(["method", "effect"])
    summary.to_csv(args.summary_output, sep="\t", index=False)


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("--bundle", type=Path, required=True)
    plan.add_argument("--min-subjects", type=int, default=4)
    plan.add_argument("--subject-folds", type=Path)
    plan.add_argument("--subject-fold", type=int)
    plan.add_argument(
        "--effect", action="append", choices=("cell_type", "condition")
    )
    plan.add_argument("--output", type=Path, required=True)
    plan.set_defaults(function=command_plan)

    split = commands.add_parser("split-subjects")
    split.add_argument("--bundle", type=Path, required=True)
    split.add_argument("--folds", type=int, default=2)
    split.add_argument("--seed", type=int, default=20260805)
    split.add_argument("--output", type=Path, required=True)
    split.set_defaults(function=command_split_subjects)

    scquint = commands.add_parser("scquint")
    scquint.add_argument("--bundle", type=Path, required=True)
    scquint.add_argument("--contrasts", type=Path, required=True)
    scquint.add_argument("--contrast-index", type=int, required=True)
    scquint.add_argument("--scquint-root", type=Path)
    scquint.add_argument("--min-samples", type=int, default=3)
    scquint.add_argument("--min-global-proportion", type=float, default=1e-3)
    scquint.add_argument("--jobs", type=int, default=1)
    scquint.add_argument("--output", type=Path, required=True)
    scquint.set_defaults(function=command_scquint)

    metadata = commands.add_parser("leafcutter-metadata")
    metadata.add_argument("--bundle", type=Path, required=True)
    metadata.add_argument("--contrasts", type=Path, required=True)
    metadata.add_argument("--contrast-index", type=int, required=True)
    metadata.add_argument("--output", type=Path, required=True)
    metadata.set_defaults(function=command_leafcutter_metadata)

    normalize = commands.add_parser("normalize-leafcutter")
    normalize.add_argument("--contrasts", type=Path, required=True)
    normalize.add_argument("--contrast-index", type=int, required=True)
    normalize.add_argument("--input", type=Path, required=True)
    normalize.add_argument("--output", type=Path, required=True)
    normalize.set_defaults(function=command_normalize_leafcutter)

    normalize_majiq = commands.add_parser("normalize-majiq")
    normalize_majiq.add_argument("--contrasts", type=Path, required=True)
    normalize_majiq.add_argument("--contrast-index", type=int, required=True)
    normalize_majiq.add_argument("--input", type=Path, required=True)
    normalize_majiq.add_argument("--output", type=Path, required=True)
    normalize_majiq.add_argument(
        "--majiq-pvalue-mode",
        choices=("raw", "posterior_q95"),
        default="raw",
    )
    normalize_majiq.set_defaults(function=command_normalize_majiq)

    normalize_majiq_dir = commands.add_parser("normalize-majiq-directory")
    normalize_majiq_dir.add_argument("--contrasts", type=Path, required=True)
    normalize_majiq_dir.add_argument("--input-dir", type=Path, required=True)
    normalize_majiq_dir.add_argument("--output-dir", type=Path, required=True)
    normalize_majiq_dir.add_argument(
        "--majiq-pvalue-mode",
        choices=("raw", "posterior_q95"),
        default="raw",
    )
    normalize_majiq_dir.set_defaults(function=command_normalize_majiq_directory)

    summarize = commands.add_parser("summarize")
    summarize.add_argument("--input", action="append", type=Path)
    summarize.add_argument("--input-dir", action="append", type=Path)
    summarize.add_argument("--output", type=Path, required=True)
    summarize.add_argument("--omnibus-output", type=Path, required=True)
    summarize.add_argument("--summary-output", type=Path, required=True)
    summarize.set_defaults(function=command_summarize)
    return result


def main():
    args = parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
