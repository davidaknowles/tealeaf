#!/usr/bin/env python3
"""Compare significant fitted local-path changes with read-derived evidence."""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, pearsonr, spearmanr

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, read_primer_cell_groups, summarize_subject_feature_evidence


def parse_interval(value):
    if pd.isna(value) or str(value).lower() == "null":
        return None
    return tuple(ast.literal_eval(str(value)))


def load_path_usage(root):
    paths = sorted(Path(root).glob("shard_*/path_usage.tsv"))
    if not paths:
        raise FileNotFoundError(f"no path-usage shards under {root}")
    usage = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    keys = ["test_id", "subject", "cell_type", "path_number"]
    if usage.duplicated(keys).any():
        raise ValueError("duplicate subject, cell-type, and path estimates")
    return usage


def build_event(record):
    anchors = tuple(interval for interval in (parse_interval(record.left_anchor), parse_interval(record.right_anchor)) if interval is not None)
    local_paths = ast.literal_eval(record.tested_path_signatures)
    paths = tuple(tuple(sorted({*anchors, *(tuple(exon) for exon in path)})) for path in local_paths)
    intervals = [interval for path in paths for interval in path]
    if not anchors or not intervals:
        return None
    splice_sites = set()
    for path in paths:
        ordered = sorted(path, reverse=record.strand == "-")
        for left, right in zip(ordered, ordered[1:]):
            junction = (left[1], right[0]) if record.strand == "+" else (right[1], left[0])
            if junction[0] < junction[1]:
                splice_sites.update(junction)
    event = SashimiEvent(record.test_id, record.chromosome, min(start for start, _ in intervals), max(end for _, end in intervals), record.strand, frozenset(splice_sites))
    return event, paths, anchors


def merge_support_by_event(results):
    exon_by_event, junction_by_event = defaultdict(Counter), defaultdict(Counter)
    for exon_blocks, junctions, _ in results:
        for key, count in exon_blocks.items():
            exon_by_event[key[0]][key] += count
        for key, count in junctions.items():
            junction_by_event[key[0]][key] += count
    return exon_by_event, junction_by_event


def make_feature_contrasts(evidence, min_common_subjects=3, min_anchor_depth=1.0, min_junction_count=5.0, min_fitted_delta=0.05):
    """Contrast raw evidence for the strongest verifiable fitted cell-type pair."""
    evidence = evidence.copy()
    evidence["adequate_support"] = np.where(evidence["feature_type"] == "exon", evidence["evidence_denominator"] >= min_anchor_depth, evidence["evidence_denominator"] >= min_junction_count)
    rows = []
    group_columns = ["test_id", "feature_id", "feature_type", "start", "end", "path_numbers"]
    for keys, feature in evidence.groupby(group_columns, sort=False, dropna=False):
        for primer, primer_data in feature[feature["adequate_support"] & feature["observed_inclusion"].notna()].groupby("primer"):
            candidates = []
            wide = primer_data.set_index(["subject", "cell_type"])[["fitted_inclusion", "observed_inclusion", "evidence_denominator"]].unstack("cell_type")
            for level_a, level_b in combinations(sorted(primer_data["cell_type"].unique()), 2):
                paired = wide.loc[:, pd.IndexSlice[:, [level_a, level_b]]].dropna()
                if len(paired) < min_common_subjects:
                    continue
                fitted_difference = float(np.mean(paired["fitted_inclusion", level_b] - paired["fitted_inclusion", level_a]))
                if fitted_difference < 0:
                    level_a, level_b = level_b, level_a
                    fitted_difference = -fitted_difference
                candidates.append((fitted_difference, level_a, level_b, paired))
            if not candidates:
                continue
            fitted_delta, low_cell_type, high_cell_type, paired = max(candidates, key=lambda value: (value[0], value[1], value[2]))
            if fitted_delta < min_fitted_delta:
                continue
            observed_differences = paired["observed_inclusion", high_cell_type].to_numpy() - paired["observed_inclusion", low_cell_type].to_numpy()
            minimum_denominator = min(paired["evidence_denominator", low_cell_type].sum(), paired["evidence_denominator", high_cell_type].sum())
            rows.append(dict(zip(group_columns, keys)) | {"primer": primer, "low_cell_type": low_cell_type, "high_cell_type": high_cell_type, "n_common_subjects": len(paired), "fitted_delta": fitted_delta, "observed_delta": float(np.mean(observed_differences)), "observed_delta_se": float(np.std(observed_differences, ddof=1) / np.sqrt(len(paired))), "minimum_total_denominator": float(minimum_denominator)})
    return pd.DataFrame(rows)


def summarize_block_contrasts(contrasts):
    return contrasts.groupby(["test_id", "feature_type", "primer"], as_index=False).agg(n_features=("feature_id", "nunique"), fitted_delta=("fitted_delta", "mean"), observed_delta=("observed_delta", "mean"))


def metric_summary(data):
    if len(data) < 3:
        return {"pearson": np.nan, "spearman": np.nan, "sign_agreement": np.nan}
    variable = data["fitted_delta"].nunique() > 1 and data["observed_delta"].nunique() > 1
    return {"pearson": pearsonr(data["fitted_delta"], data["observed_delta"]).statistic if variable else np.nan, "spearman": spearmanr(data["fitted_delta"], data["observed_delta"]).statistic if variable else np.nan, "sign_agreement": float((data["observed_delta"] > 0).mean())}


def summarize_agreement(blocks, bootstrap_replicates=1000, seed=13):
    rng = np.random.default_rng(seed)
    rows = []
    for (feature_type, primer), group in blocks.groupby(["feature_type", "primer"], sort=True):
        observed = metric_summary(group)
        bootstrap = defaultdict(list)
        for _ in range(bootstrap_replicates):
            sample = group.iloc[rng.integers(0, len(group), len(group))]
            for metric, value in metric_summary(sample).items():
                bootstrap[metric].append(value)
        row = {"feature_type": feature_type, "primer": primer, "n_blocks": group["test_id"].nunique(), "n_features": int(group["n_features"].sum()), **observed}
        for metric, values in bootstrap.items():
            row[f"{metric}_ci_low"], row[f"{metric}_ci_high"] = np.nanquantile(values, [0.025, 0.975])
        rows.append(row)
    return pd.DataFrame(rows)


def write_concordance_tables(evidence, output_dir):
    output_dir = Path(output_dir)
    contrasts = make_feature_contrasts(evidence)
    blocks = summarize_block_contrasts(contrasts)
    summary = summarize_agreement(blocks)
    contrasts.to_csv(output_dir / "feature_contrasts.tsv", sep="\t", index=False, na_rep="NA")
    blocks.to_csv(output_dir / "block_contrasts.tsv", sep="\t", index=False, na_rep="NA")
    summary.to_csv(output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    rank_features, rank_blocks, rank_summary = summarize_rank_concordance(evidence)
    rank_features.to_csv(output_dir / "feature_rank_concordance.tsv", sep="\t", index=False, na_rep="NA")
    rank_blocks.to_csv(output_dir / "block_rank_concordance.tsv", sep="\t", index=False, na_rep="NA")
    rank_summary.to_csv(output_dir / "rank_summary.tsv", sep="\t", index=False, na_rep="NA")
    return contrasts, blocks, summary


def summarize_rank_concordance(evidence, min_total_denominator=20.0, min_cell_types=3, min_fitted_range=0.05, bootstrap_replicates=1000, seed=17):
    """Rank fitted and raw cell-type means within features, then blocks."""
    data = evidence.copy()
    if "endpoint_depth" not in data:
        group_columns = ["test_id", "subject", "cell_type", "primer"]
        endpoints = data[data["feature_type"] == "exon"].groupby(group_columns)["evidence_denominator"].first().rename("endpoint_depth")
        data = data.join(endpoints, on=group_columns)
    if "feature_signal" not in data:
        data["feature_signal"] = np.where(data["feature_type"] == "exon", data["unclipped_observed_ratio"] * data["evidence_denominator"], data["observed_inclusion"] * data["evidence_denominator"])
    group_columns = ["test_id", "feature_id", "feature_type", "primer", "cell_type"]
    pooled = data.groupby(group_columns, as_index=False).agg(fitted_inclusion=("fitted_inclusion", "mean"), feature_signal=("feature_signal", "sum"), endpoint_depth=("endpoint_depth", "sum"))
    pooled["observed_inclusion"] = pooled["feature_signal"] / pooled["endpoint_depth"]
    pooled = pooled[(pooled["endpoint_depth"] >= min_total_denominator) & pooled["observed_inclusion"].notna()]
    rows = []
    for keys, feature in pooled.groupby(["test_id", "feature_id", "feature_type", "primer"], sort=False):
        fitted_range = feature["fitted_inclusion"].max() - feature["fitted_inclusion"].min()
        if len(feature) < min_cell_types or fitted_range < min_fitted_range or feature["observed_inclusion"].nunique() < 2:
            continue
        rows.append(dict(zip(["test_id", "feature_id", "feature_type", "primer"], keys)) | {"n_cell_types": len(feature), "fitted_range": fitted_range, "rank_correlation": spearmanr(feature["fitted_inclusion"], feature["observed_inclusion"]).statistic})
    features = pd.DataFrame(rows)
    blocks = features.groupby(["test_id", "feature_type", "primer"], as_index=False).agg(n_features=("feature_id", "nunique"), rank_correlation=("rank_correlation", "mean"))
    rng = np.random.default_rng(seed)
    summaries = []
    for (feature_type, primer), group in blocks.groupby(["feature_type", "primer"], sort=True):
        bootstrap_median, bootstrap_positive = [], []
        for _ in range(bootstrap_replicates):
            sample = group.iloc[rng.integers(0, len(group), len(group))]["rank_correlation"]
            bootstrap_median.append(sample.median())
            bootstrap_positive.append((sample > 0).mean())
        positive = int((group["rank_correlation"] > 0).sum())
        summaries.append({"feature_type": feature_type, "primer": primer, "n_blocks": len(group), "n_features": int(group["n_features"].sum()), "median_rank_correlation": group["rank_correlation"].median(), "median_ci_low": np.quantile(bootstrap_median, 0.025), "median_ci_high": np.quantile(bootstrap_median, 0.975), "positive_fraction": positive / len(group), "positive_ci_low": np.quantile(bootstrap_positive, 0.025), "positive_ci_high": np.quantile(bootstrap_positive, 0.975), "positive_sign_p_value": binomtest(positive, len(group), 0.5, alternative="greater").pvalue})
    return features, blocks, pd.DataFrame(summaries)


def write_plot(output_dir):
    from plotnine import aes, coord_cartesian, element_blank, element_text, facet_grid, geom_abline, geom_boxplot, geom_hline, geom_jitter, geom_point, geom_smooth, geom_text, ggplot, labs, scale_color_manual, theme, theme_bw

    output_dir = Path(output_dir)
    blocks = pd.read_csv(output_dir / "block_contrasts.tsv", sep="\t")
    summary = pd.read_csv(output_dir / "summary.tsv", sep="\t")
    label_map = {"exon": "Variable exon segments", "junction": "Variable junctions"}
    blocks["feature_class"] = blocks["feature_type"].map(label_map)
    blocks["primer"] = pd.Categorical(blocks["primer"], categories=["poly(dT)", "random hexamer"], ordered=True)
    annotations = summary.copy()
    annotations["feature_class"] = annotations["feature_type"].map(label_map)
    annotations["label"] = annotations.apply(lambda row: f"{int(row.n_blocks)} blocks\nSpearman {row.spearman:.2f} [{row.spearman_ci_low:.2f}, {row.spearman_ci_high:.2f}]\nraw direction agrees {row.sign_agreement:.0%} [{row.sign_agreement_ci_low:.0%}, {row.sign_agreement_ci_high:.0%}]", axis=1)
    annotations["x"] = 0.03
    annotations["y"] = 0.94
    plot = ggplot(blocks, aes("fitted_delta", "observed_delta", color="feature_class"))
    plot += geom_hline(yintercept=0, color="#777777", size=0.45)
    plot += geom_abline(intercept=0, slope=1, linetype="dashed", color="#777777", size=0.45)
    plot += geom_smooth(method="lm", se=True, size=0.75)
    plot += geom_point(alpha=0.58, size=1.35)
    plot += geom_text(data=annotations, mapping=aes("x", "y", label="label"), inherit_aes=False, ha="left", va="top", size=7)
    plot += facet_grid("feature_class ~ primer")
    plot += scale_color_manual(values={"Variable exon segments": "#3B6FB6", "Variable junctions": "#C5523C"})
    plot += coord_cartesian(xlim=(0, 1), ylim=(-1, 1))
    plot += labs(x="Inferred high-minus-low feature inclusion", y="Read-derived high-minus-low inclusion", title="Endpoint-normalized read evidence for significant local-path changes", caption="Each point averages the variable features of one BH-significant block. A cell-type pair is selected from fitted usage among pairs with common subjects, without inspecting raw evidence. Dashed line, equality; ribbon, linear fit with 95% confidence interval.")
    plot += theme_bw(base_size=10)
    plot += theme(legend_position="none", panel_grid_minor=element_blank(), strip_text=element_text(size=9), plot_title=element_text(size=11), plot_caption=element_text(size=8))
    plot.save(output_dir / "path_coverage_concordance.pdf", width=8.3, height=6.4, units="in", verbose=False)

    ranks = pd.read_csv(output_dir / "block_rank_concordance.tsv", sep="\t")
    rank_summary = pd.read_csv(output_dir / "rank_summary.tsv", sep="\t")
    ranks["feature_class"] = ranks["feature_type"].map(label_map)
    ranks["primer"] = pd.Categorical(ranks["primer"], categories=["poly(dT)", "random hexamer"], ordered=True)
    rank_summary["feature_class"] = rank_summary["feature_type"].map(label_map)
    rank_summary["primer"] = pd.Categorical(rank_summary["primer"], categories=["poly(dT)", "random hexamer"], ordered=True)
    rank_summary["label"] = rank_summary.apply(lambda row: f"{int(row.n_blocks)} blocks\nmedian {row.median_rank_correlation:.2f} [{row.median_ci_low:.2f}, {row.median_ci_high:.2f}]\npositive {row.positive_fraction:.0%} [{row.positive_ci_low:.0%}, {row.positive_ci_high:.0%}]", axis=1)
    rank_summary["y"] = 1.19
    rank_plot = ggplot(ranks, aes("primer", "rank_correlation", color="primer"))
    rank_plot += geom_hline(yintercept=0, color="#555555", size=0.55)
    rank_plot += geom_boxplot(width=0.14, outlier_alpha=0, fill="white", size=0.65)
    rank_plot += geom_jitter(width=0.2, height=0, alpha=1.0, size=0.55)
    rank_plot += geom_text(data=rank_summary, mapping=aes("primer", "y", label="label"), inherit_aes=False, va="top", size=7)
    rank_plot += facet_grid(". ~ feature_class")
    rank_plot += scale_color_manual(values={"poly(dT)": "#0B6666", "random hexamer": "#A94E16"})
    rank_plot += coord_cartesian(ylim=(-1, 1.25))
    rank_plot += labs(x=None, y="Within-block rank agreement", title="Fitted path usage and endpoint-normalized read evidence rank cell types similarly")
    rank_plot += theme_bw(base_size=10)
    rank_plot += theme(legend_position="none", panel_grid_minor=element_blank(), strip_text=element_text(size=9), plot_title=element_text(size=11), plot_caption=element_text(size=8))
    rank_plot.save(output_dir / "path_coverage_rank_concordance.pdf", width=8.3, height=4.8, units="in", verbose=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starsolo-root", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--event-catalog", type=Path)
    parser.add_argument("--path-usage-root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--bam-name", default="spliced_pseudobulk.bam")
    parser.add_argument("--deduplicate", action="store_true")
    parser.add_argument("--skip-plot", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    if args.plot_only:
        write_plot(args.output_dir)
        return
    if args.summarize_only:
        evidence = pd.read_csv(args.output_dir / "subject_feature_evidence.tsv.gz", sep="\t")
        write_concordance_tables(evidence, args.output_dir)
        if not args.skip_plot:
            write_plot(args.output_dir)
        return
    required = (args.starsolo_root, args.metadata, args.primer_pairs, args.event_catalog, args.path_usage_root)
    if any(value is None for value in required):
        parser.error("all data inputs are required unless --plot-only is used")
    catalog = pd.read_csv(args.event_catalog, sep="\t")
    catalog = catalog[pd.to_numeric(catalog["fdr"], errors="coerce") < 0.05].copy()
    usage = load_path_usage(args.path_usage_root)
    usage = usage[usage["test_id"].isin(catalog["test_id"])].copy()
    event_data = {}
    for record in catalog.itertuples(index=False):
        built = build_event(record)
        if built is not None:
            event_data[record.test_id] = built
    events = tuple(value[0] for value in event_data.values())
    cell_types = tuple(sorted(usage["cell_type"].unique()))
    barcode_groups, _ = read_primer_cell_groups(args.metadata, args.primer_pairs, cell_types)
    bam_paths = sorted(args.starsolo_root.glob(f"*/{args.bam_name}"))
    if not bam_paths:
        raise FileNotFoundError(f"no spliced pseudobulk BAMs under {args.starsolo_root}")
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(collect_bam_event_support, path, barcode_groups, events, args.deduplicate) for path in bam_paths]
        exon_by_event, junction_by_event = merge_support_by_event(future.result() for future in as_completed(futures))
    evidence_tables = []
    for test_id, (event, paths, anchors) in event_data.items():
        local_usage = usage[usage["test_id"] == test_id]
        evidence = summarize_subject_feature_evidence(event, paths, anchors, exon_by_event[test_id], junction_by_event[test_id], local_usage)
        if not evidence.empty:
            evidence_tables.append(evidence)
    evidence = pd.concat(evidence_tables, ignore_index=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(args.output_dir / "subject_feature_evidence.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    write_concordance_tables(evidence, args.output_dir)
    manifest = {"selection": "all BH-significant full-data omnibus local blocks with at least one explicit endpoint anchor", "n_significant_blocks": int(len(catalog)), "n_blocks_with_endpoint_anchors": int(len(event_data)), "n_bam_files": len(bam_paths), "bam_name": args.bam_name, "deduplicated_in_collector": args.deduplicate, "reporting_concentration": 1.0, "endpoint_normalized_evidence": "variable atomic-exon aligned depth or variable-junction count divided by mean aligned depth over the explicit block endpoint anchor exon or exons", "splice_site_psi": "the complementary contrast table uses junction count divided by the total count over unique observed junctions sharing its donor or acceptor splice site", "contrast": "for each variable feature and primer, choose the largest fitted difference among cell-type pairs with at least three common support-qualified subjects, orient it high-minus-low without inspecting raw evidence, and then calculate the raw contrast", "rank_concordance": "rank cell types by fitted feature inclusion and endpoint-normalized feature signal; require at least three cell types, total endpoint depth at least 20, and fitted range at least 0.05; average feature Spearman correlations within each block", "support": "endpoint depth at least one for each exon observation, junction denominator at least five for each junction observation, and paired fitted difference at least 0.05", "interpretation": "same-read descriptive concordance check, not independent validation or an input to statistical calibration"}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not args.skip_plot:
        write_plot(args.output_dir)


if __name__ == "__main__":
    main()
