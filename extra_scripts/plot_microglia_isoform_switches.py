#!/usr/bin/env python3
"""Find fitted dominant-path switches and plot primer-separated PSI evidence."""

from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, combine_path_usage_and_coverage, merge_support, read_primer_cell_groups, summarize_path_ordered_psi, summarize_path_usage, write_path_evidence_heatmap


def _parse_interval(value):
    if pd.isna(value) or str(value).lower() == "null":
        return None
    return tuple(ast.literal_eval(str(value)))


def load_path_usage(root):
    """Load subject-level path estimates from all completed fit shards."""
    paths = sorted(Path(root).glob("shard_*/path_usage.tsv"))
    if not paths:
        raise FileNotFoundError(f"no shard path-usage tables under {root}")
    usage = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    if usage.duplicated(["test_id", "subject", "cell_type", "path_number"]).any():
        raise ValueError("duplicate subject, cell-type, and path estimates across shards")
    return usage


def rank_switches(path_usage, catalog, fdr_threshold=0.05):
    """Rank significant blocks by the largest cell-type mean path-usage span."""
    means = summarize_path_usage(path_usage)
    dominant = means.loc[means.groupby(["test_id", "cell_type"])["mean_proportion"].idxmax(), ["test_id", "cell_type", "path_number"]]
    dominant_counts = dominant.groupby("test_id")["path_number"].nunique().rename("n_dominant_paths")
    extrema = means.groupby(["test_id", "path_number", "path"], as_index=False)["mean_proportion"].agg(min_usage="min", max_usage="max")
    extrema["usage_span"] = extrema["max_usage"] - extrema["min_usage"]
    extrema = extrema.sort_values(["test_id", "usage_span", "path_number"], ascending=[True, False, True]).drop_duplicates("test_id")
    min_rows = means.loc[means.groupby(["test_id", "path_number"])["mean_proportion"].idxmin(), ["test_id", "path_number", "cell_type"]].rename(columns={"cell_type": "min_cell_type"})
    max_rows = means.loc[means.groupby(["test_id", "path_number"])["mean_proportion"].idxmax(), ["test_id", "path_number", "cell_type"]].rename(columns={"cell_type": "max_cell_type"})
    annotations = catalog.loc[pd.to_numeric(catalog["fdr"], errors="coerce") < fdr_threshold].copy()
    ranking = extrema.merge(min_rows, on=["test_id", "path_number"]).merge(max_rows, on=["test_id", "path_number"]).merge(dominant_counts, on="test_id").merge(annotations, on="test_id", validate="one_to_one")
    ranking["dominant_path_switch"] = ranking["n_dominant_paths"] > 1
    return ranking.sort_values(["dominant_path_switch", "usage_span", "test_id"], ascending=[False, False, True]).reset_index(drop=True), means


def _event_from_row(row):
    anchors = [interval for interval in (_parse_interval(row.left_anchor), _parse_interval(row.right_anchor)) if interval is not None]
    local_paths = ast.literal_eval(row.tested_path_signatures)
    paths = [sorted({*anchors, *(tuple(exon) for exon in local_path)}) for local_path in local_paths]
    intervals = [interval for path in paths for interval in path]
    if not intervals:
        raise ValueError(f"no intervals for {row.test_id}")
    splice_sites = set()
    for path in paths:
        ordered = sorted(path, reverse=row.strand == "-")
        for left, right in zip(ordered, ordered[1:]):
            junction = (left[1], right[0]) if row.strand == "+" else (right[1], left[0])
            if junction[0] < junction[1]:
                splice_sites.update(junction)
    event = SashimiEvent(row.gene_name, row.chromosome, min(start for start, _ in intervals), max(end for _, end in intervals), row.strand, frozenset(splice_sites))
    return event, paths


def _has_spliced_junction(row):
    event, paths = _event_from_row(row)
    for path in paths:
        ordered = sorted(path, reverse=event.strand == "-")
        for left, right in zip(ordered, ordered[1:]):
            junction = (left[1], right[0]) if event.strand == "+" else (right[1], left[0])
            if junction[0] < junction[1]:
                return True
    return False


def render_selected(output_dir):
    """Render saved switch matrices without reopening BAM files."""
    selected = pd.read_csv(output_dir / "selected_switches.tsv", sep="\t")
    for row in selected.itertuples(index=False):
        matrix = pd.read_csv(output_dir / row.gene_name / f"{row.gene_name}_psi_heatmap.tsv", sep="\t")
        cell_types = tuple(matrix["cell_type"].drop_duplicates())
        for primer in ("poly(dT)", "random hexamer"):
            write_path_evidence_heatmap(matrix, row.gene_name, primer, cell_types, output_dir, evidence_name="Exon or junction PSI", evidence_limits=(0, 1), output_suffix="psi_heatmap")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starsolo-root", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--event-catalog", type=Path)
    parser.add_argument("--path-usage-root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--n-events", type=int, default=3)
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    if args.plot_only:
        render_selected(args.output_dir)
        return
    required = (args.starsolo_root, args.metadata, args.primer_pairs, args.event_catalog, args.path_usage_root)
    if any(value is None for value in required):
        parser.error("BAM, metadata, catalog, and path-usage inputs are required unless --plot-only is used")
    usage = load_path_usage(args.path_usage_root)
    catalog = pd.read_csv(args.event_catalog, sep="\t")
    ranking, usage_summary = rank_switches(usage, catalog)
    ranking["has_spliced_junction"] = ranking.apply(_has_spliced_junction, axis=1)
    interpretable = ranking[ranking["dominant_path_switch"] & ranking["has_spliced_junction"] & ~ranking["event_type"].str.contains("whole-gene")]
    selected = interpretable.head(args.n_events).copy()
    if len(selected) < args.n_events:
        raise ValueError(f"only {len(selected)} significant dominant-path switches were found")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(args.output_dir / "switch_ranking.tsv", sep="\t", index=False)
    selected.to_csv(args.output_dir / "selected_switches.tsv", sep="\t", index=False)
    selected_usage = usage[usage["test_id"].isin(selected["test_id"])]
    selected_usage.to_csv(args.output_dir / "path_usage.tsv", sep="\t", index=False)
    usage_summary[usage_summary["test_id"].isin(selected["test_id"])].to_csv(args.output_dir / "path_usage_summary.tsv", sep="\t", index=False)
    events, paths_by_gene = [], {}
    for row in selected.itertuples(index=False):
        event, paths = _event_from_row(row)
        events.append(event)
        paths_by_gene[event.event_id] = paths
    cell_types = tuple(sorted(selected_usage["cell_type"].unique()))
    barcode_groups, cell_sizes = read_primer_cell_groups(args.metadata, args.primer_pairs, cell_types)
    bam_paths = sorted(args.starsolo_root.glob("*/spliced_pseudobulk.bam"))
    if not bam_paths:
        raise FileNotFoundError(f"no spliced pseudobulk BAMs under {args.starsolo_root}")
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(collect_bam_event_support, bam_paths, repeat(barcode_groups), repeat(tuple(events))))
    exon_blocks, junctions, event_umis = merge_support(results)
    cell_sizes.to_csv(args.output_dir / "cell_counts.tsv", sep="\t", index=False)
    pd.DataFrame([{"event": key[0], "subject": key[1], "cell_type": key[2], "primer": key[3], "spliced_umis": value} for key, value in sorted(event_umis.items())]).to_csv(args.output_dir / "spliced_umi_support.tsv", sep="\t", index=False)
    for row in selected.itertuples(index=False):
        event = next(event for event in events if event.event_id == row.gene_name)
        event_cell_types = tuple(sorted(selected_usage.loc[selected_usage["test_id"] == row.test_id, "cell_type"].unique()))
        psi = summarize_path_ordered_psi(event, paths_by_gene[event.event_id], exon_blocks, junctions, event_cell_types)
        matrix = combine_path_usage_and_coverage(usage_summary, psi, row.test_id)
        event_dir = args.output_dir / event.event_id
        event_dir.mkdir(parents=True, exist_ok=True)
        matrix.to_csv(event_dir / f"{event.event_id}_psi_heatmap.tsv", sep="\t", index=False)
    manifest = {
        "selection": "three non-whole-gene BH-significant blocks with the largest fitted cell-type mean path-usage span among blocks whose dominant fitted path changes between cell types and whose displayed paths contain at least one splice junction",
        "literal_zero_to_one_switches": int(((ranking["min_usage"] <= 0.1) & (ranking["max_usage"] >= 0.9)).sum()),
        "exon_psi": "mean aligned depth divided by the mean depth of the nearest path-shared exon on each available side, clipped to [0,1]; path-shared exons are 1",
        "junction_psi": "junction UMI count divided by the total UMI count of unique observed junctions sharing its donor or acceptor splice site",
        "events": selected[["gene_name", "test_id", "event_type", "path_number", "min_cell_type", "min_usage", "max_cell_type", "max_usage", "usage_span", "fdr"]].to_dict("records"),
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not args.skip_plots:
        render_selected(args.output_dir)


if __name__ == "__main__":
    main()
