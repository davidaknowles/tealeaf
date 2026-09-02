#!/usr/bin/env python3
"""Prepare priming-separated read support for selected microglia events."""

from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import json
from pathlib import Path
import subprocess

import pandas as pd

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, combine_path_usage_and_coverage, merge_support, read_primer_cell_groups, select_strongest_significant_contrasts, summarize_path_ordered_coverage, summarize_path_usage, write_ggsashimi_inputs


EVENTS = (
    ("App", "ENSMUSG00000022892.12:B3"),
    ("Gria2", "ENSMUSG00000033981.15:B1"),
    ("Grin1", "ENSMUSG00000026959.14:B3"),
)


def _event_gtf(path, event, row):
    paths = ast.literal_eval(row.tested_path_signatures)
    anchors = [ast.literal_eval(row.left_anchor), ast.literal_eval(row.right_anchor)]
    lines = []
    for index, internal_exons in enumerate(paths, start=1):
        exons = sorted({tuple(exon) for exon in [*anchors, *internal_exons]})
        transcript_id = f"{event.event_id}.path{index}"
        attributes = f'gene_id "{event.event_id}"; transcript_id "{transcript_id}";'
        lines.append("\t".join((event.chromosome, "Tealeaf", "transcript", str(min(start for start, _ in exons) + 1), str(max(end for _, end in exons)), ".", event.strand, ".", attributes)))
        for exon_number, (start, end) in enumerate(exons, start=1):
            exon_attributes = attributes + f' exon_number "{exon_number}";'
            lines.append("\t".join((event.chromosome, "Tealeaf", "exon", str(start + 1), str(end), ".", event.strand, ".", exon_attributes)))
    path.write_text("\n".join(lines) + "\n")


def write_path_usage_plot(summary, event_id, test_id, cell_types, output_dir):
    """Plot mean fitted path usage with subject-level standard errors."""
    from plotnine import aes, element_text, geom_col, geom_errorbar, ggplot, labs, position_dodge, scale_fill_manual, scale_x_discrete, scale_y_continuous, theme, theme_bw

    event_summary = summary[summary["test_id"] == test_id].copy()
    if event_summary.empty:
        raise ValueError(f"no fitted path usage for {test_id}")
    event_summary["cell_type"] = pd.Categorical(event_summary["cell_type"], categories=cell_types, ordered=True)
    dodge = position_dodge(width=0.78)
    plot = (
        ggplot(event_summary, aes("cell_type", "mean_proportion", fill="path"))
        + geom_col(position=dodge, width=0.72, color="white")
        + geom_errorbar(aes(ymin="mean_proportion - se_proportion", ymax="mean_proportion + se_proportion"), position=dodge, width=0.15)
        + scale_fill_manual(values=["#7570B3", "#E7298A"])
        + scale_x_discrete(labels=lambda values: [value.replace("_", " ") for value in values])
        + scale_y_continuous(limits=(0, 1), breaks=[0, 0.25, 0.5, 0.75, 1], labels=lambda values: [f"{value:.0%}" for value in values], expand=(0, 0.02))
        + labs(x=None, y="Estimated path usage", fill=None)
        + theme_bw(base_size=11)
        + theme(axis_text_x=element_text(rotation=15, ha="right"), legend_position="top", figure_size=(7.2, 2.8))
    )
    output_path = output_dir / event_id / f"{event_id}_path_usage.pdf"
    plot.save(output_path, width=7.2, height=2.8, units="in", verbose=False)
    return output_path


def write_block_heatmap(matrix, event_id, primer, cell_types, output_dir):
    """Plot one primer-specific, path-ordered block evidence heatmap."""
    from plotnine import aes, element_blank, element_text, geom_text, geom_tile, geom_vline, ggplot, labs, scale_fill_cmap, scale_x_discrete, theme, theme_bw

    data = matrix[matrix["primer"] == primer].copy()
    if data.empty:
        raise ValueError(f"no {primer} evidence for {event_id}")
    columns = data[["feature_id", "feature", "feature_type", "column_order"]].drop_duplicates().sort_values("column_order")
    feature_ids = columns["feature_id"].tolist()
    feature_labels = dict(zip(columns["feature_id"], columns["feature"]))
    data["feature_id"] = pd.Categorical(data["feature_id"], categories=feature_ids, ordered=True)
    data["cell_type"] = pd.Categorical(data["cell_type"], categories=tuple(cell_types)[::-1], ordered=True)
    data["coverage"] = data["raw_value"].where(data["feature_type"] != "path")
    data["label"] = data.apply(lambda row: f"{row.raw_value:.0%}" if row.feature_type == "path" else "", axis=1)
    path_positions = [index for index, value in enumerate(columns["feature_type"], start=1) if value == "path"]
    primer_slug = "oligodt" if primer == "poly(dT)" else "ranhex"
    plot = (
        ggplot(data, aes("feature_id", "cell_type"))
        + geom_tile(aes(fill="coverage"), color="white", size=0.15)
        + geom_text(aes(label="label"), size=6)
        + geom_vline(xintercept=[position - 0.5 for position in path_positions[1:]], color="#25364A", size=0.7)
        + scale_fill_cmap(name="Coverage per 1,000 half-cells", cmap_name="viridis", na_value="#F2F2F2")
        + scale_x_discrete(labels=lambda values: [feature_labels[value] for value in values])
        + labs(x="Path and ordered block components", y=None, title=f"{event_id}, {primer}")
        + theme_bw(base_size=9)
        + theme(axis_text_x=element_text(size=7, rotation=90, va="top"), axis_text_y=element_text(size=7), panel_grid=element_blank(), legend_position="right", plot_title=element_text(size=10))
    )
    output_path = output_dir / event_id / f"{event_id}_{primer_slug}_block_heatmap.pdf"
    height = max(3.0, 1.4 + 0.32 * len(cell_types))
    plot.save(output_path, width=9.5, height=height, units="in", verbose=False)
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starsolo-root", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--primer-pairs", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--pairwise-event-catalog", type=Path, required=True)
    parser.add_argument("--pairwise-fit-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--ggsashimi", type=Path)
    parser.add_argument("--run-plots", action="store_true")
    parser.add_argument("--skip-heatmap-plots", action="store_true")
    parser.add_argument("--skip-path-usage-plot", action="store_true")
    parser.add_argument("--path-usage", type=Path)
    parser.add_argument("--omnibus-path-usage", type=Path, required=True)
    args = parser.parse_args()
    omnibus_path_usage = pd.read_csv(args.omnibus_path_usage, sep="\t")
    omnibus_path_usage_summary = summarize_path_usage(omnibus_path_usage)
    heatmap_levels = {block_id: tuple(sorted(group["cell_type"].unique())) for block_id, group in omnibus_path_usage.groupby("block_id")}
    block_ids = [event_id for _, event_id in EVENTS]
    selected_contrasts = select_strongest_significant_contrasts(sorted(args.pairwise_fit_root.glob("*/paired_path.tsv")), args.pairwise_event_catalog, block_ids)
    contrasts = {row.block_id: (row.level_a, row.level_b) for row in selected_contrasts.itertuples(index=False)}
    catalog = pd.read_csv(args.event_catalog, sep="\t", dtype=str)
    rows, events = {}, []
    for gene, event_id in EVENTS:
        cell_types = contrasts[event_id]
        matched = catalog[(catalog["gene_name"] == gene) & (catalog["test_id"] == event_id)]
        if len(matched) != 1:
            raise ValueError(f"expected one catalog row for {gene} {event_id}, found {len(matched)}")
        row = matched.iloc[0]
        anchors = [ast.literal_eval(row.left_anchor), ast.literal_eval(row.right_anchor)]
        paths = ast.literal_eval(row.tested_path_signatures)
        intervals = [*anchors, *(exon for local_path in paths for exon in local_path)]
        events.append(SashimiEvent(gene, row.chromosome, min(interval[0] for interval in intervals), max(interval[1] for interval in intervals), row.strand))
        rows[gene] = row
        rows[gene + "_cell_types"] = cell_types
        rows[gene + "_paths"] = [sorted({tuple(exon) for exon in [*anchors, *local_path]}) for local_path in paths]
        rows[gene + "_heatmap_cell_types"] = heatmap_levels[event_id]
    selected_cell_types = {cell_type for levels in heatmap_levels.values() for cell_type in levels}
    barcode_groups, cell_sizes = read_primer_cell_groups(args.metadata, args.primer_pairs, selected_cell_types)
    bam_paths = sorted(args.starsolo_root.glob("*/spliced_pseudobulk.bam"))
    if not bam_paths:
        raise FileNotFoundError(f"no spliced pseudobulk BAMs under {args.starsolo_root}")
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(collect_bam_event_support, bam_paths, repeat(barcode_groups), repeat(tuple(events))))
    exon_blocks, junctions, event_umis = merge_support(results)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    omnibus_path_usage.to_csv(args.output_dir / "omnibus_path_usage.tsv", sep="\t", index=False)
    omnibus_path_usage_summary.to_csv(args.output_dir / "omnibus_path_usage_summary.tsv", sep="\t", index=False)
    path_usage_summary = None
    if args.path_usage is not None:
        path_usage = pd.read_csv(args.path_usage, sep="\t")
        path_usage_summary = summarize_path_usage(path_usage)
        path_usage.to_csv(args.output_dir / "path_usage.tsv", sep="\t", index=False)
        path_usage_summary.to_csv(args.output_dir / "path_usage_summary.tsv", sep="\t", index=False)
    cell_sizes.to_csv(args.output_dir / "cell_counts.tsv", sep="\t", index=False)
    selected_contrasts.to_csv(args.output_dir / "selected_contrasts.tsv", sep="\t", index=False)
    support_rows = []
    for key, count in sorted(event_umis.items()):
        event_id, subject, cell_type, primer = key
        support_rows.append({"event": event_id, "subject": subject, "cell_type": cell_type, "primer": primer, "spliced_umis": count})
    pd.DataFrame(support_rows).to_csv(args.output_dir / "spliced_umi_support.tsv", sep="\t", index=False)
    manifest = {"contrast_selection": "largest mean_difference_norm among pairwise contrasts with BH FDR < 0.05", "normalization": "UMI-deduplicated spliced alignments per 1000 primer-specific half-cells", "block_heatmaps": "Path markers followed by each path's exon and junction components in 5-prime-to-3-prime RNA order; shared components are repeated under each path; coverage colors are raw primer-specific values per 1000 half-cells; path columns show fitted percentages", "runs": len(bam_paths), "events": []}
    palette = args.output_dir / "palette.tsv"
    palette.write_text("#1B9E77\n#D95F02\n#1B9E77\n#D95F02\n")
    event_ids = dict(EVENTS)
    for event in events:
        event_dir = args.output_dir / event.event_id
        cell_types = rows[event.event_id + "_cell_types"]
        intron, exon, strand = write_ggsashimi_inputs(event_dir, event, exon_blocks, junctions, cell_sizes, cell_types)
        gtf = event_dir / f"{event.event_id}_paths.gtf"
        _event_gtf(gtf, event, rows[event.event_id])
        heatmap_cell_types = rows[event.event_id + "_heatmap_cell_types"]
        coverage = summarize_path_ordered_coverage(event, rows[event.event_id + "_paths"], exon_blocks, junctions, cell_sizes, heatmap_cell_types)
        matrix = combine_path_usage_and_coverage(omnibus_path_usage_summary, coverage, event_ids[event.event_id])
        matrix.to_csv(event_dir / f"{event.event_id}_block_heatmap.tsv", sep="\t", index=False)
        oligodt_heatmap = event_dir / f"{event.event_id}_oligodt_block_heatmap.pdf"
        ranhex_heatmap = event_dir / f"{event.event_id}_ranhex_block_heatmap.pdf"
        if not args.skip_heatmap_plots:
            write_block_heatmap(matrix, event.event_id, "poly(dT)", heatmap_cell_types, args.output_dir)
            write_block_heatmap(matrix, event.event_id, "random hexamer", heatmap_cell_types, args.output_dir)
        coordinates = f"{event.chromosome}:{event.start + 1}-{event.end + 1}"
        prefix = event_dir / f"{event.event_id}_priming_sashimi"
        contrast = selected_contrasts[selected_contrasts["block_id"] == event_ids[event.event_id]].iloc[0]
        manifest["events"].append({"gene": event.event_id, "test_id": event_ids[event.event_id], "cell_types": list(cell_types), "mean_difference_norm": float(contrast.mean_difference_norm), "pairwise_fdr": float(contrast.fdr), "coordinates": coordinates, "intron_counts": str(intron.relative_to(args.output_dir)), "exon_counts": str(exon.relative_to(args.output_dir)), "path_annotation": str(gtf.relative_to(args.output_dir)), "figure": str(prefix.with_suffix('.pdf').relative_to(args.output_dir))})
        manifest["events"][-1]["heatmap_cell_types"] = list(heatmap_cell_types)
        manifest["events"][-1]["block_heatmap_data"] = str((event_dir / f"{event.event_id}_block_heatmap.tsv").relative_to(args.output_dir))
        if path_usage_summary is not None:
            usage_figure = event_dir / f"{event.event_id}_path_usage.pdf"
            if not args.skip_path_usage_plot:
                write_path_usage_plot(path_usage_summary, event.event_id, contrast.test_id, cell_types, args.output_dir)
            manifest["events"][-1]["path_usage_figure"] = str(usage_figure.relative_to(args.output_dir))
        manifest["events"][-1]["oligodt_block_heatmap"] = str(oligodt_heatmap.relative_to(args.output_dir))
        manifest["events"][-1]["ranhex_block_heatmap"] = str(ranhex_heatmap.relative_to(args.output_dir))
        if args.run_plots:
            if args.ggsashimi is None:
                raise ValueError("--ggsashimi is required with --run-plots")
            subprocess.run(["python", str(args.ggsashimi), "--intron", str(intron), "--exon", str(exon), "--strand_info", str(strand), "--data_type", "sc", "--aggregation", "False", "-c", coordinates, "-g", str(gtf), "-o", str(prefix), "-M", "1", "--shrink", "--fix-y-scale", "--height", "1.1", "--ann-height", "0.9", "--width", "9", "--base-size", "11", "-P", str(palette), "-F", "pdf"], check=True)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
