"""Read-level support summaries for single-cell sashimi plots."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pysam
except ModuleNotFoundError:
    pysam = None

from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode


@dataclass(frozen=True)
class SashimiEvent:
    """A zero-based, half-open genomic interval to summarize."""

    event_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    splice_sites: frozenset[int] | None = None


def read_primer_cell_groups(metadata_path, primer_pairs_path, cell_types):
    """Map half-cell barcodes to subject, cell type, and primer."""
    primer_pairs = pd.read_csv(primer_pairs_path, sep="\t", dtype=str)
    primer_by_barcode = {}
    for primer, column in (("poly(dT)", "polydt_barcode"), ("random hexamer", "ranhex_barcode")):
        for barcode in primer_pairs[column].dropna():
            if barcode in primer_by_barcode:
                raise ValueError(f"barcode {barcode!r} occurs in both primer groups")
            primer_by_barcode[barcode] = primer
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    metadata = metadata[metadata["cell_type"].isin(cell_types)].copy()
    metadata["primer"] = metadata["barcode"].map(primer_by_barcode)
    metadata = metadata.dropna(subset=["primer"])
    if metadata.empty:
        raise ValueError("no selected cells have a primer assignment")
    if metadata["barcode"].duplicated().any():
        raise ValueError("selected metadata contains duplicate barcodes")
    lookup = {row.barcode: (row.subject, row.cell_type, row.primer) for row in metadata.itertuples(index=False)}
    sizes = metadata.groupby(["subject", "cell_type", "primer"], sort=True).size().rename("cells").reset_index()
    return lookup, sizes


def select_strongest_significant_contrasts(fit_paths, event_catalog_path, block_ids, fdr_threshold=0.05):
    """Select each block's largest fitted effect among BH-significant contrasts."""
    columns = ["test_id", "block_id", "level_a", "level_b", "mean_difference_norm"]
    fit_tables = [pd.read_csv(path, sep="\t", usecols=columns) for path in fit_paths]
    if not fit_tables:
        raise ValueError("no pairwise fit tables were provided")
    fits = pd.concat(fit_tables, ignore_index=True)
    if fits["test_id"].duplicated().any():
        raise ValueError("pairwise fit tables contain duplicate test identifiers")
    catalog_columns = ["test_id", "block_id", "statistic", "empirical_p_value", "fdr"]
    catalog = pd.read_csv(event_catalog_path, sep="\t", usecols=catalog_columns)
    catalog["fdr"] = pd.to_numeric(catalog["fdr"], errors="coerce")
    significant = catalog[(catalog["block_id"].isin(block_ids)) & (catalog["fdr"] < fdr_threshold)]
    candidates = significant.merge(fits, on=["test_id", "block_id"], how="inner", validate="one_to_one")
    missing_fits = sorted(set(significant["test_id"]) - set(candidates["test_id"]))
    if missing_fits:
        raise ValueError(f"significant contrasts without fitted effects: {missing_fits}")
    candidates["mean_difference_norm"] = pd.to_numeric(candidates["mean_difference_norm"], errors="coerce")
    candidates = candidates.dropna(subset=["mean_difference_norm"])
    selected = candidates.sort_values(["block_id", "mean_difference_norm", "test_id"], ascending=[True, False, True]).drop_duplicates("block_id")
    missing = sorted(set(block_ids) - set(selected["block_id"]))
    if missing:
        raise ValueError(f"blocks without a significant fitted contrast: {missing}")
    return selected.sort_values("block_id").reset_index(drop=True)


def summarize_path_usage(path_usage):
    """Summarize subject-level fitted path proportions as means and standard errors."""
    required = {"test_id", "cell_type", "path", "path_number", "proportion"}
    missing = required - set(path_usage.columns)
    if missing:
        raise ValueError(f"path-usage table is missing columns: {sorted(missing)}")
    summary = path_usage.groupby(["test_id", "cell_type", "path", "path_number"], sort=True)["proportion"].agg(n_subjects="count", mean_proportion="mean", sd_proportion="std").reset_index()
    summary["sd_proportion"] = summary["sd_proportion"].fillna(0.0)
    summary["se_proportion"] = summary["sd_proportion"] / np.sqrt(summary["n_subjects"])
    return summary


def variable_path_features(event, paths):
    """Return variable exon segments and junctions with their path membership.

    Exons are split at every boundary before path membership is assigned, so an
    alternative donor or acceptor contributes only its discriminating segment.
    """
    ordered_paths = [tuple(sorted(set(tuple(exon) for exon in path), reverse=event.strand == "-")) for path in paths]
    n_paths = len(ordered_paths)
    if n_paths < 2:
        raise ValueError("at least two paths are required")
    boundaries = sorted({coordinate for path in ordered_paths for exon in path for coordinate in exon})
    exon_segments = []
    for start, end in zip(boundaries, boundaries[1:]):
        membership = tuple(index + 1 for index, path in enumerate(ordered_paths) if any(exon_start <= start and end <= exon_end for exon_start, exon_end in path))
        if membership and len(membership) < n_paths:
            if exon_segments and exon_segments[-1][1] == start and exon_segments[-1][2] == membership:
                exon_segments[-1] = (exon_segments[-1][0], end, membership)
            else:
                exon_segments.append((start, end, membership))

    def junction_between(left_exon, right_exon):
        return (left_exon[1], right_exon[0]) if event.strand == "+" else (right_exon[1], left_exon[0])

    junction_membership = defaultdict(list)
    for path_number, path in enumerate(ordered_paths, start=1):
        for left, right in zip(path, path[1:]):
            start, end = junction_between(left, right)
            if start < end:
                junction_membership[(start, end)].append(path_number)
    rows = []
    for start, end, membership in exon_segments:
        rows.append({"feature_type": "exon", "start": start, "end": end, "path_numbers": membership})
    for (start, end), membership in sorted(junction_membership.items()):
        if len(membership) < n_paths:
            rows.append({"feature_type": "junction", "start": start, "end": end, "path_numbers": tuple(membership)})
    features = pd.DataFrame(rows)
    if features.empty:
        return pd.DataFrame(columns=["feature_type", "start", "end", "path_numbers", "feature_id"])
    features = features.sort_values(["feature_type", "start", "end"]).reset_index(drop=True)
    features["feature_id"] = features.apply(lambda row: f"{row.feature_type[0].upper()}:{int(row.start)}-{int(row.end)}", axis=1)
    return features


def summarize_subject_feature_evidence(event, paths, anchors, exon_blocks, junctions, path_usage, primers=("poly(dT)", "random hexamer")):
    """Compare fitted path inclusion with subject-level read evidence.

    For each subject, cell type, primer, and variable feature, the fitted value
    is the sum of the proportions of paths containing that feature. Exon depth
    is divided by mean depth across the explicit block endpoint anchors and
    clipped to [0, 1]. Junction counts are divided by counts over unique
    junctions sharing their donor or acceptor site.
    """
    required = {"subject", "cell_type", "path_number", "proportion"}
    missing = required - set(path_usage.columns)
    if missing:
        raise ValueError(f"path usage is missing columns: {sorted(missing)}")
    anchors = tuple(tuple(anchor) for anchor in anchors)
    if not anchors:
        raise ValueError("at least one explicit endpoint anchor is required")
    features = variable_path_features(event, paths)
    if features.empty:
        return pd.DataFrame()
    exon_by_group = defaultdict(list)
    for key, count in exon_blocks.items():
        event_id, subject, cell_type, primer, start, end = key
        if event_id == event.event_id:
            exon_by_group[(subject, cell_type, primer)].append((start, end, count))
    junction_by_group = defaultdict(dict)
    for key, count in junctions.items():
        event_id, subject, cell_type, primer, start, end = key
        if event_id == event.event_id:
            group_counts = junction_by_group[(subject, cell_type, primer)]
            group_counts[(start, end)] = group_counts.get((start, end), 0.0) + count

    def mean_depth(blocks, interval):
        start, end = interval
        aligned_bases = sum(max(0, min(block_end, end) - max(block_start, start)) * count for block_start, block_end, count in blocks)
        return aligned_bases / (end - start)

    rows = []
    usage_groups = path_usage.groupby(["subject", "cell_type"], sort=True)
    for (subject, cell_type), group_usage in usage_groups:
        proportions = dict(zip(group_usage["path_number"].astype(int), group_usage["proportion"].astype(float)))
        for primer in primers:
            group = (subject, cell_type, primer)
            blocks = exon_by_group.get(group, ())
            counts = junction_by_group.get(group, {})
            anchor_depth = float(np.mean([mean_depth(blocks, anchor) for anchor in anchors]))
            for feature in features.itertuples(index=False):
                fitted = sum(proportions.get(path_number, 0.0) for path_number in feature.path_numbers)
                if feature.feature_type == "exon":
                    feature_signal = mean_depth(blocks, (feature.start, feature.end))
                    raw_ratio = feature_signal / anchor_depth if anchor_depth > 0 else np.nan
                    observed = np.clip(raw_ratio, 0.0, 1.0) if np.isfinite(raw_ratio) else np.nan
                    denominator = anchor_depth
                else:
                    feature_signal = counts.get((feature.start, feature.end), 0.0)
                    denominator = sum(count for (start, end), count in counts.items() if start == feature.start or end == feature.end)
                    raw_ratio = feature_signal / denominator if denominator > 0 else np.nan
                    observed = raw_ratio
                endpoint_ratio = feature_signal / anchor_depth if anchor_depth > 0 else np.nan
                rows.append({"test_id": event.event_id, "subject": subject, "cell_type": cell_type, "primer": primer, "feature_id": feature.feature_id, "feature_type": feature.feature_type, "start": feature.start, "end": feature.end, "path_numbers": ",".join(map(str, feature.path_numbers)), "fitted_inclusion": fitted, "observed_inclusion": observed, "unclipped_observed_ratio": raw_ratio, "evidence_denominator": denominator, "feature_signal": feature_signal, "endpoint_depth": anchor_depth, "endpoint_normalized_signal": endpoint_ratio})
    return pd.DataFrame(rows)


def coverage_group_columns(cell_types):
    """Map cell-type and primer labels to ggsashimi coverage columns."""
    primers = ("poly(dT)", "random hexamer")

    def compact(value):
        return value.translate(str.maketrans("", "", "_ ()"))
    return {(cell_type, primer): f"{compact(cell_type)}_{compact(primer)}_0" for cell_type in cell_types for primer in primers}


def summarize_path_ordered_coverage(event, paths, exon_blocks, junctions, cell_sizes, cell_types, scale=1000.0):
    """Summarize path-ordered exon and junction coverage by cell type and primer."""
    groups = coverage_group_columns(cell_types)
    denominators = cell_sizes.groupby(["cell_type", "primer"])["cells"].sum().to_dict()
    missing = [group for group in groups if denominators.get(group, 0) == 0]
    if missing:
        raise ValueError(f"groups without cells: {missing}")
    ordered_paths = [sorted(set(tuple(exon) for exon in path), reverse=event.strand == "-") for path in paths]
    exon_order = sorted({exon for path in ordered_paths for exon in path}, reverse=event.strand == "-")
    exon_labels = {exon: f"E{index}" for index, exon in enumerate(exon_order, start=1)}

    def junction_between(left_exon, right_exon):
        left_start, left_end = left_exon
        right_start, right_end = right_exon
        return (left_end, right_start) if event.strand == "+" else (right_end, left_start)

    junction_order = sorted({junction_between(left, right) for path in ordered_paths for left, right in zip(path, path[1:])}, reverse=event.strand == "-")
    junction_labels = {junction: f"J{index}" for index, junction in enumerate(junction_order, start=1)}
    rows = []
    column_order = 0
    for path_number, ordered_exons in enumerate(ordered_paths, start=1):
        column_order += 1
        for exon_index, (exon_start, exon_end) in enumerate(ordered_exons):
            feature = exon_labels[(exon_start, exon_end)]
            feature_id = f"P{path_number}:{feature}"
            coordinate = f"{event.chromosome}:{exon_start + 1}-{exon_end}"
            for cell_type, primer in groups:
                aligned_bases = 0.0
                for key, count in exon_blocks.items():
                    event_id, _, local_cell_type, local_primer, start, end = key
                    if event_id == event.event_id and local_cell_type == cell_type and local_primer == primer:
                        aligned_bases += max(0, min(end, exon_end) - max(start, exon_start)) * count
                coverage = aligned_bases * scale / (denominators[(cell_type, primer)] * (exon_end - exon_start))
                rows.append({"path_number": path_number, "feature_type": "exon", "feature": feature, "feature_id": feature_id, "coordinate": coordinate, "column_order": column_order, "cell_type": cell_type, "primer": primer, "raw_value": coverage})
            column_order += 1
            if exon_index == len(ordered_exons) - 1:
                continue
            junction_start, junction_end = junction_between(ordered_exons[exon_index], ordered_exons[exon_index + 1])
            feature = junction_labels[(junction_start, junction_end)]
            feature_id = f"P{path_number}:{feature}"
            coordinate = f"{event.chromosome}:{junction_start + 1}-{junction_end + 1}"
            for cell_type, primer in groups:
                count = sum(value for key, value in junctions.items() if key[0] == event.event_id and key[2] == cell_type and key[3] == primer and key[4] == junction_start and key[5] == junction_end)
                coverage = count * scale / denominators[(cell_type, primer)]
                rows.append({"path_number": path_number, "feature_type": "junction", "feature": feature, "feature_id": feature_id, "coordinate": coordinate, "column_order": column_order, "cell_type": cell_type, "primer": primer, "raw_value": coverage})
            column_order += 1
    return pd.DataFrame(rows)


def summarize_path_ordered_psi(event, paths, exon_blocks, junctions, cell_types):
    """Summarize exon and splice-site-normalized junction PSI along tested paths."""
    groups = tuple(coverage_group_columns(cell_types))
    ordered_paths = [sorted(set(tuple(exon) for exon in path), reverse=event.strand == "-") for path in paths]
    exon_order = sorted({exon for path in ordered_paths for exon in path}, reverse=event.strand == "-")
    exon_labels = {exon: f"E{index}" for index, exon in enumerate(exon_order, start=1)}
    shared_exons = set.intersection(*(set(path) for path in ordered_paths))

    def junction_between(left_exon, right_exon):
        left_start, left_end = left_exon
        right_start, right_end = right_exon
        return (left_end, right_start) if event.strand == "+" else (right_end, left_start)

    path_junctions = [{junction_between(left, right) for left, right in zip(path, path[1:]) if junction_between(left, right)[0] < junction_between(left, right)[1]} for path in ordered_paths]
    shared_junctions = set.intersection(*path_junctions)
    junction_order = sorted(set.union(*path_junctions), reverse=event.strand == "-")
    junction_labels = {junction: f"J{index}" for index, junction in enumerate(junction_order, start=1)}
    exon_depths = {}
    junction_counts = {}
    for cell_type, primer in groups:
        for exon_start, exon_end in exon_order:
            aligned_bases = sum(max(0, min(end, exon_end) - max(start, exon_start)) * count for key, count in exon_blocks.items() for event_id, _, local_cell_type, local_primer, start, end in (key,) if event_id == event.event_id and local_cell_type == cell_type and local_primer == primer)
            exon_depths[(cell_type, primer, (exon_start, exon_end))] = aligned_bases / (exon_end - exon_start)
        counts = {}
        for key, count in junctions.items():
            event_id, _, local_cell_type, local_primer, start, end = key
            if event_id == event.event_id and local_cell_type == cell_type and local_primer == primer:
                counts[(start, end)] = counts.get((start, end), 0.0) + count
        junction_counts[(cell_type, primer)] = counts
    rows = []
    column_order = 0
    exon_positions = {exon: index for index, exon in enumerate(exon_order)}
    shared_positions = sorted(exon_positions[exon] for exon in shared_exons)
    for path_number, ordered_exons in enumerate(ordered_paths, start=1):
        column_order += 1
        for cell_type, primer in groups:
            rows.append({"path_number": path_number, "feature_type": "path_order", "feature": f"P{path_number}", "feature_id": f"P{path_number}:order", "coordinate": "", "column_order": column_order - 1, "cell_type": cell_type, "primer": primer, "raw_value": np.nan})
        for exon_index, exon in enumerate(ordered_exons):
            if exon not in shared_exons:
                feature = exon_labels[exon]
                coordinate = f"{event.chromosome}:{exon[0] + 1}-{exon[1]}"
                position = exon_positions[exon]
                left = [index for index in shared_positions if index < position]
                right = [index for index in shared_positions if index > position]
                reference_exons = [exon_order[index] for index in ([max(left)] if left else []) + ([min(right)] if right else [])]
                for cell_type, primer in groups:
                    depth = exon_depths[(cell_type, primer, exon)]
                    reference_depths = [exon_depths[(cell_type, primer, reference)] for reference in reference_exons]
                    reference_depth = float(np.mean(reference_depths)) if reference_depths else np.nan
                    psi = np.clip(depth / reference_depth, 0.0, 1.0) if reference_depth > 0 else np.nan
                    rows.append({"path_number": path_number, "feature_type": "exon", "feature": feature, "feature_id": f"P{path_number}:{feature}", "coordinate": coordinate, "column_order": column_order, "cell_type": cell_type, "primer": primer, "raw_value": psi})
                column_order += 1
            if exon_index == len(ordered_exons) - 1:
                continue
            junction = junction_between(ordered_exons[exon_index], ordered_exons[exon_index + 1])
            if junction[0] >= junction[1] or junction in shared_junctions:
                continue
            feature = junction_labels[junction]
            coordinate = f"{event.chromosome}:{junction[0] + 1}-{junction[1] + 1}"
            for cell_type, primer in groups:
                counts = junction_counts[(cell_type, primer)]
                denominator = sum(count for alternative, count in counts.items() if alternative[0] == junction[0] or alternative[1] == junction[1])
                psi = counts.get(junction, 0.0) / denominator if denominator > 0 else np.nan
                rows.append({"path_number": path_number, "feature_type": "junction", "feature": feature, "feature_id": f"P{path_number}:{feature}", "coordinate": coordinate, "column_order": column_order, "cell_type": cell_type, "primer": primer, "raw_value": psi})
            column_order += 1
    return pd.DataFrame(rows)


def combine_path_usage_and_coverage(path_usage_summary, coverage, test_id):
    """Interleave fitted path usage with path-ordered coverage features."""
    usage = path_usage_summary[path_usage_summary["test_id"] == test_id]
    if usage.empty:
        raise ValueError(f"no fitted path usage for {test_id}")
    markers = []
    primers = tuple(coverage["primer"].drop_duplicates())
    path_orders = coverage[coverage["feature_type"] == "path_order"]
    evidence = coverage[coverage["feature_type"] != "path_order"]
    for path_number in sorted(usage["path_number"].unique()):
        path_usage = usage[usage["path_number"] == path_number]
        explicit_order = path_orders[path_orders["path_number"] == path_number]
        if not explicit_order.empty:
            marker_order = int(explicit_order["column_order"].iloc[0])
        else:
            marker_order = int(evidence.loc[evidence["path_number"] == path_number, "column_order"].min()) - 1
        for primer in primers:
            for record in path_usage.itertuples(index=False):
                markers.append({"path_number": path_number, "feature_type": "path", "feature": f"P{path_number}", "feature_id": f"P{path_number}:path", "coordinate": "", "column_order": marker_order, "cell_type": record.cell_type, "primer": primer, "raw_value": float(record.mean_proportion)})
    combined = pd.concat([pd.DataFrame(markers), evidence], ignore_index=True)
    return combined.sort_values(["primer", "column_order", "cell_type"]).reset_index(drop=True)


def write_path_evidence_heatmap(matrix, event_id, primer, cell_types, output_dir, evidence_name="Coverage per 1,000 half-cells", evidence_limits=None, output_suffix="block_heatmap", shared_scale=False):
    """Plot one primer-specific matrix of path usage and ordered evidence."""
    from plotnine import aes, element_blank, element_text, geom_point, geom_text, geom_tile, geom_vline, ggplot, labs, scale_color_cmap, scale_fill_cmap, scale_x_discrete, theme, theme_bw

    data = matrix[matrix["primer"] == primer].copy()
    if data.empty:
        raise ValueError(f"no {primer} evidence for {event_id}")
    columns = data[["feature_id", "feature", "feature_type", "column_order"]].drop_duplicates().sort_values("column_order")
    feature_ids = columns["feature_id"].tolist()
    feature_labels = dict(zip(columns["feature_id"], columns["feature"]))
    data["feature_id"] = pd.Categorical(data["feature_id"], categories=feature_ids, ordered=True)
    data["cell_type"] = pd.Categorical(data["cell_type"], categories=tuple(cell_types)[::-1], ordered=True)
    data["evidence"] = data["raw_value"].where(data["feature_type"] != "path")
    data["label"] = data.apply(lambda row: f"{row.raw_value:.0%}" if row.feature_type == "path" else "", axis=1)
    path_data = data[data["feature_type"] == "path"]
    path_positions = [index for index, value in enumerate(columns["feature_type"], start=1) if value == "path"]
    marker_size = min(14, 270 / len(feature_ids))
    primer_slug = "oligodt" if primer == "poly(dT)" else "ranhex"
    fill_column = "raw_value" if shared_scale else "evidence"
    plot = ggplot(data, aes("feature_id", "cell_type")) + geom_tile(aes(fill=fill_column), color="white", size=0.15)
    if not shared_scale:
        plot += geom_point(data=path_data, mapping=aes(color="raw_value"), shape="s", size=marker_size)
        plot += scale_color_cmap(name="Path usage", cmap_name="magma", limits=(0, 1))
    plot += geom_text(data=path_data, mapping=aes(label="label"), size=6, fontweight="bold")
    plot += geom_vline(xintercept=[position - 0.5 for position in path_positions[1:]], color="#25364A", size=0.7)
    plot += scale_fill_cmap(name=evidence_name, cmap_name="viridis", limits=evidence_limits, na_value="#F2F2F2")
    plot += scale_x_discrete(labels=lambda values: [feature_labels[value] for value in values])
    plot += labs(x="Path and ordered block components", y=None, title=f"{event_id}, {primer}")
    plot += theme_bw(base_size=9)
    plot += theme(axis_text_x=element_text(size=7, rotation=90, va="top"), axis_text_y=element_text(size=7), panel_grid=element_blank(), legend_position="right", plot_title=element_text(size=10))
    output_path = Path(output_dir) / event_id / f"{event_id}_{primer_slug}_{output_suffix}.pdf"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    height = max(3.0, 1.4 + 0.32 * len(cell_types))
    plot.save(output_path, width=9.5, height=height, units="in", verbose=False)
    return output_path


def _alignment_junctions(alignment):
    position = alignment.reference_start
    for operation, length in alignment.cigartuples or ():
        if operation == 3:
            yield position, position + length
        if operation in {0, 2, 3, 7, 8}:
            position += length


def collect_bam_event_support(bam_path, barcode_groups, events, deduplicate=False):
    """Collect aligned blocks and junctions, optionally deduplicating cell UMIs."""
    if pysam is None:
        raise ImportError("pysam is required to collect BAM event support")
    bin_size = 1_000_000
    by_chromosome = defaultdict(lambda: defaultdict(list))
    for event in events:
        for bin_index in range(event.start // bin_size, (event.end - 1) // bin_size + 1):
            by_chromosome[event.chromosome][bin_index].append(event)
    exon_blocks, junctions, event_umis = Counter(), Counter(), Counter()
    seen = set()
    with pysam.AlignmentFile(bam_path, "rb") as alignments:
        for alignment in alignments.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
            chromosome_bins = by_chromosome.get(alignment.reference_name)
            if chromosome_bins is None or not alignment.has_tag("CB") or alignment.reference_end is None:
                continue
            local_events = {}
            for bin_index in range(alignment.reference_start // bin_size, (alignment.reference_end - 1) // bin_size + 1):
                local_events.update((event.event_id, event) for event in chromosome_bins.get(bin_index, ()))
            overlapping = [event for event in local_events.values() if alignment.reference_start < event.end and alignment.reference_end > event.start]
            if not overlapping:
                continue
            barcode = normalize_starsolo_barcode(alignment.get_tag("CB"))
            group = barcode_groups.get(barcode)
            if group is None:
                continue
            if deduplicate:
                if not alignment.has_tag("NH") or alignment.get_tag("NH") != 1 or not alignment.has_tag("UB"):
                    continue
                key = (barcode, alignment.get_tag("UB"), alignment.reference_id, alignment.reference_start, alignment.cigarstring)
                if key in seen:
                    continue
                seen.add(key)
            subject, cell_type, primer = group
            for event in overlapping:
                event_umis[(event.event_id, subject, cell_type, primer)] += 1
                for start, end in alignment.get_blocks():
                    start, end = max(start, event.start), min(end, event.end)
                    if start < end:
                        exon_blocks[(event.event_id, subject, cell_type, primer, start, end)] += 1
                for start, end in _alignment_junctions(alignment):
                    within_interval = event.start <= start and end <= event.end
                    shares_target_site = event.splice_sites is not None and (start in event.splice_sites or end in event.splice_sites)
                    if shares_target_site or (event.splice_sites is None and within_interval):
                        junctions[(event.event_id, subject, cell_type, primer, start, end)] += 1
    return exon_blocks, junctions, event_umis


def merge_support(results):
    """Add support counters returned for independent BAMs."""
    merged = Counter(), Counter(), Counter()
    for result in results:
        for target, source in zip(merged, result):
            target.update(source)
    return merged


def _constant_segments(values):
    changes = np.flatnonzero(np.diff(values) != 0) + 1
    boundaries = np.concatenate(([0], changes, [len(values)]))
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        value = float(values[start])
        if value > 0:
            yield int(start), int(end), value


def write_ggsashimi_inputs(output_dir, event, exon_blocks, junctions, cell_sizes, cell_types, scale=1000.0):
    """Write ggsashimi tables normalized per scale primer-specific half-cells."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    group_names = coverage_group_columns(cell_types)
    denominators = cell_sizes.groupby(["cell_type", "primer"])["cells"].sum().to_dict()
    missing = [group for group in group_names if denominators.get(group, 0) == 0]
    if missing:
        raise ValueError(f"groups without cells: {missing}")
    exon_rows, junction_rows = {}, {}
    for group, name in group_names.items():
        cell_type, primer = group
        denominator = denominators[group]
        coverage = np.zeros(event.end - event.start, dtype=float)
        for key, count in exon_blocks.items():
            event_id, _, local_cell_type, local_primer, start, end = key
            if event_id == event.event_id and local_cell_type == cell_type and local_primer == primer:
                coverage[start - event.start : end - event.start] += count
        coverage *= scale / denominator
        for start, end, value in _constant_segments(coverage):
            genomic_start, genomic_end = event.start + start + 1, event.start + end + 1
            row = exon_rows.setdefault((genomic_start, genomic_end), {"Name": f"{event.chromosome}:{genomic_start}-{genomic_end}", "Chr": event.chromosome, "Start": genomic_start, "End": genomic_end})
            row[name] = value
        for key, count in junctions.items():
            event_id, _, local_cell_type, local_primer, start, end = key
            if event_id == event.event_id and local_cell_type == cell_type and local_primer == primer:
                genomic_start, genomic_end = start + 1, end + 1
                row = junction_rows.setdefault((genomic_start, genomic_end), {"Name": f"{event.chromosome}:{genomic_start}-{genomic_end}", "Chr": event.chromosome, "Start": genomic_start, "End": genomic_end})
                row[name] = row.get(name, 0.0) + count * scale / denominator
    columns = ["Name", "Chr", "Start", "End", *group_names.values()]
    for rows in (exon_rows, junction_rows):
        for row in rows.values():
            for column in group_names.values():
                row.setdefault(column, 0.0)
    exon_path = output_dir / f"{event.event_id}_exon.tsv"
    intron_path = output_dir / f"{event.event_id}_intron.tsv"
    strand_path = output_dir / f"{event.event_id}_strand.tsv"
    pd.DataFrame(exon_rows.values(), columns=columns).sort_values(["Start", "End"]).to_csv(exon_path, sep=" ", index=False)
    pd.DataFrame(junction_rows.values(), columns=columns).sort_values(["Start", "End"]).to_csv(intron_path, sep=" ", index=False)
    with open(strand_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter=" ", lineterminator="\n")
        writer.writerow(("intron", "near_exons", "strand"))
        writer.writerow(("dummy", "dummy", event.strand))
    return intron_path, exon_path, strand_path
