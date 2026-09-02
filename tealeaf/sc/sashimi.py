"""Read-level support summaries for single-cell sashimi plots."""

from __future__ import annotations

from collections import Counter
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

    junction_order = sorted({junction_between(left, right) for path in ordered_paths for left, right in zip(path, path[1:]) if junction_between(left, right)[0] < junction_between(left, right)[1]}, reverse=event.strand == "-")
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
        for exon_index, exon in enumerate(ordered_exons):
            feature = exon_labels[exon]
            coordinate = f"{event.chromosome}:{exon[0] + 1}-{exon[1]}"
            position = exon_positions[exon]
            left = [index for index in shared_positions if index < position]
            right = [index for index in shared_positions if index > position]
            reference_exons = [exon_order[index] for index in ([max(left)] if left else []) + ([min(right)] if right else [])]
            for cell_type, primer in groups:
                depth = exon_depths[(cell_type, primer, exon)]
                if exon in shared_exons:
                    psi = 1.0
                else:
                    reference_depths = [exon_depths[(cell_type, primer, reference)] for reference in reference_exons]
                    reference_depth = float(np.mean(reference_depths)) if reference_depths else np.nan
                    psi = np.clip(depth / reference_depth, 0.0, 1.0) if reference_depth > 0 else np.nan
                rows.append({"path_number": path_number, "feature_type": "exon", "feature": feature, "feature_id": f"P{path_number}:{feature}", "coordinate": coordinate, "column_order": column_order, "cell_type": cell_type, "primer": primer, "raw_value": psi})
            column_order += 1
            if exon_index == len(ordered_exons) - 1:
                continue
            junction = junction_between(ordered_exons[exon_index], ordered_exons[exon_index + 1])
            if junction[0] >= junction[1]:
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
    for path_number in sorted(coverage["path_number"].unique()):
        path_usage = usage[usage["path_number"] == path_number]
        marker_order = int(coverage.loc[coverage["path_number"] == path_number, "column_order"].min()) - 1
        for primer in primers:
            for record in path_usage.itertuples(index=False):
                markers.append({"path_number": path_number, "feature_type": "path", "feature": f"P{path_number}", "feature_id": f"P{path_number}:path", "coordinate": "", "column_order": marker_order, "cell_type": record.cell_type, "primer": primer, "raw_value": float(record.mean_proportion)})
    combined = pd.concat([pd.DataFrame(markers), coverage], ignore_index=True)
    return combined.sort_values(["primer", "column_order", "cell_type"]).reset_index(drop=True)


def write_path_evidence_heatmap(matrix, event_id, primer, cell_types, output_dir, evidence_name="Coverage per 1,000 half-cells", evidence_limits=None, output_suffix="block_heatmap"):
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
    plot = (
        ggplot(data, aes("feature_id", "cell_type"))
        + geom_tile(aes(fill="evidence"), color="white", size=0.15)
        + geom_point(data=path_data, mapping=aes(color="raw_value"), shape="s", size=marker_size)
        + geom_text(data=path_data, mapping=aes(label="label"), size=6, fontweight="bold")
        + geom_vline(xintercept=[position - 0.5 for position in path_positions[1:]], color="#25364A", size=0.7)
        + scale_fill_cmap(name=evidence_name, cmap_name="viridis", limits=evidence_limits, na_value="#F2F2F2")
        + scale_color_cmap(name="Path usage", cmap_name="magma", limits=(0, 1))
        + scale_x_discrete(labels=lambda values: [feature_labels[value] for value in values])
        + labs(x="Path and ordered block components", y=None, title=f"{event_id}, {primer}")
        + theme_bw(base_size=9)
        + theme(axis_text_x=element_text(size=7, rotation=90, va="top"), axis_text_y=element_text(size=7), panel_grid=element_blank(), legend_position="right", plot_title=element_text(size=10))
    )
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


def collect_bam_event_support(bam_path, barcode_groups, events):
    """Collect UMI-deduplicated aligned blocks and junctions from one BAM."""
    if pysam is None:
        raise ImportError("pysam is required to collect BAM event support")
    by_chromosome = {}
    for event in events:
        by_chromosome.setdefault(event.chromosome, []).append(event)
    exon_blocks, junctions, event_umis = Counter(), Counter(), Counter()
    with pysam.AlignmentFile(bam_path, "rb") as alignments:
        for alignment in alignments.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
            local_events = by_chromosome.get(alignment.reference_name, ())
            if not local_events or not alignment.has_tag("CB") or alignment.reference_end is None:
                continue
            overlapping = [event for event in local_events if alignment.reference_start < event.end and alignment.reference_end > event.start]
            if not overlapping:
                continue
            barcode = normalize_starsolo_barcode(alignment.get_tag("CB"))
            group = barcode_groups.get(barcode)
            if group is None:
                continue
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
