"""Read-level support summaries for single-cell sashimi plots."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import csv
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode


@dataclass(frozen=True)
class SashimiEvent:
    """A zero-based, half-open genomic interval to summarize."""

    event_id: str
    chromosome: str
    start: int
    end: int
    strand: str


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


def summarize_junction_coverage(table, group_columns):
    """Return long-form normalized coverage for ordered junction features."""
    ordered = table.sort_values(["Start", "End"]).reset_index(drop=True).copy()
    ordered["feature"] = [f"J{index}" for index in range(1, len(ordered) + 1)]
    ordered["coordinate"] = ordered["Chr"].astype(str) + ":" + ordered["Start"].astype(str) + "-" + ordered["End"].astype(str)
    rows = []
    for (cell_type, primer), column in group_columns.items():
        if column not in ordered:
            raise ValueError(f"coverage table is missing column {column!r}")
        for record in ordered[["feature", "coordinate", column]].itertuples(index=False, name=None):
            feature, coordinate, value = record
            rows.append({"feature_type": "junction", "feature": feature, "coordinate": coordinate, "cell_type": cell_type, "primer": primer, "coverage": float(value)})
    return pd.DataFrame(rows)


def summarize_exon_coverage(table, exon_intervals, group_columns):
    """Average normalized base coverage over ordered, zero-based half-open exons."""
    segments = table.copy()
    segments["segment_start"] = segments["Start"].astype(int) - 1
    segments["segment_end"] = segments["End"].astype(int) - 1
    rows = []
    for exon_index, (chromosome, exon_start, exon_end) in enumerate(sorted(set(exon_intervals), key=lambda value: (value[1], value[2])), start=1):
        feature = f"E{exon_index}"
        coordinate = f"{chromosome}:{exon_start + 1}-{exon_end}"
        overlap = np.maximum(0, np.minimum(segments["segment_end"].to_numpy(), exon_end) - np.maximum(segments["segment_start"].to_numpy(), exon_start))
        for (cell_type, primer), column in group_columns.items():
            if column not in segments:
                raise ValueError(f"coverage table is missing column {column!r}")
            coverage = float(np.dot(overlap, segments[column].to_numpy(dtype=float)) / (exon_end - exon_start))
            rows.append({"feature_type": "exon", "feature": feature, "coordinate": coordinate, "cell_type": cell_type, "primer": primer, "coverage": coverage})
    return pd.DataFrame(rows)


def _alignment_junctions(alignment):
    position = alignment.reference_start
    for operation, length in alignment.cigartuples or ():
        if operation == 3:
            yield position, position + length
        if operation in {0, 2, 3, 7, 8}:
            position += length


def collect_bam_event_support(bam_path, barcode_groups, events):
    """Collect UMI-deduplicated aligned blocks and junctions from one BAM."""
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
                    if event.start <= start and end <= event.end:
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
