#!/usr/bin/env python3
"""Score fitted log-normalized gene PCA against reference cell labels."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc import representation_scoring


def named_prefix(value):
    name, separator, prefix = value.partition("=")
    if not separator or not name or not prefix:
        raise argparse.ArgumentTypeError("fits must have the form NAME=OUTPUT_PREFIX")
    return name, Path(prefix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", action="append", required=True, type=named_prefix)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--groups", type=Path)
    parser.add_argument("--transcript-to-gene", required=True, type=Path)
    eligible = parser.add_mutually_exclusive_group()
    eligible.add_argument("--standard-h5ad", type=Path)
    eligible.add_argument(
        "--eligible-genes",
        type=Path,
        help="one gene id per line, or a TSV whose first column contains gene ids",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--pca-components", type=int, default=30)
    parser.add_argument("--hvg", type=int, default=2_000)
    parser.add_argument("--target-sum", type=float, default=10_000.0)
    parser.add_argument("--reconstruction-batch-cells", type=int, default=2_048)
    parser.add_argument("--silhouette-sample-size", type=int, default=10_000)
    parser.add_argument("--device", default="auto")
    args = parser.parse_args()

    labels = pd.read_csv(
        args.labels, header=None, names=["cell_id", "label"], dtype=str
    )
    groups = None
    if args.groups is not None:
        group_table = pd.read_csv(
            args.groups, header=None, names=["cell_id", "combined"], dtype=str
        )
        groups = group_table[["cell_id"]].assign(
            group=group_table["combined"].str.rsplit("__", n=1).str[-1]
        )

    transcript_to_gene = pd.read_csv(
        args.transcript_to_gene,
        sep="\t",
        header=None,
        names=["transcript_id", "gene_id"],
        dtype=str,
    )
    if args.standard_h5ad is not None:
        import anndata as ad

        standard = ad.read_h5ad(args.standard_h5ad, backed="r")
        if "gene_id" not in standard.var:
            raise ValueError("standard AnnData var must contain gene_id")
        eligible_genes = standard.var["gene_id"].astype(str).to_numpy()
        standard.file.close()
    elif args.eligible_genes is not None:
        eligible_genes = (
            pd.read_csv(
                args.eligible_genes,
                sep="\t",
                header=None,
                usecols=[0],
                dtype=str,
            )
            .iloc[:, 0]
            .dropna()
            .drop_duplicates()
            .to_numpy()
        )
    else:
        eligible_genes = (
            transcript_to_gene["gene_id"].dropna().drop_duplicates().to_numpy()
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for name, prefix in args.fit:
        output_prefix = args.output_dir / f"{name}_"
        left = right = None
        try:
            cell_ids, transcripts, left, right = (
                representation_scoring.load_glm_factors(prefix)
            )
            positions, aligned_labels, aligned_groups, aligned_ids = (
                representation_scoring.align_reference_metadata(
                    cell_ids, labels, groups
                )
            )
            gene_left, mapping_diagnostics = (
                representation_scoring.aggregate_transcript_loadings(
                    left, transcripts, transcript_to_gene, eligible_genes
                )
            )
            embedding, selected_genes, active, preprocessing = (
                representation_scoring.log_gene_pca_embedding(
                    right[positions],
                    gene_left,
                    eligible_genes,
                    target_sum=args.target_sum,
                    n_hvg=args.hvg,
                    n_components=args.pca_components,
                    batch_cells=args.reconstruction_batch_cells,
                    device=args.device,
                )
            )
            preprocessing.update(mapping_diagnostics)
            representation_scoring.write_embedding(
                embedding,
                aligned_ids,
                selected_genes,
                preprocessing,
                output_prefix,
            )
            report, folds = representation_scoring.score_embedding(
                embedding,
                aligned_labels,
                aligned_groups,
                valid_mask=active,
                n_splits=args.folds,
                silhouette_sample_size=args.silhouette_sample_size,
            )
            report.update(preprocessing)
            report["status"] = "ok"
        except (OSError, ValueError, RuntimeError) as error:
            report = {"status": "invalid", "error": str(error)}
            folds = pd.DataFrame()
        report.update(name=name, fit_prefix=str(prefix))
        representation_scoring.write_score(report, folds, output_prefix)
        summaries.append(report)
        print(json.dumps(report, indent=2), flush=True)
        left = right = None
        gc.collect()
    pd.DataFrame(summaries).to_csv(args.output_dir / "label_score_summary.csv", index=False)


if __name__ == "__main__":
    main()
