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

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, merge_support, read_primer_cell_groups, write_ggsashimi_inputs


EVENTS = (
    ("App", "ENSMUSG00000022892.12:B3", ("EX_cortical", "INH_midbrain")),
    ("Gria2", "ENSMUSG00000033981.15:B1", ("EX_cortical", "INH_midbrain")),
    ("Grin1", "ENSMUSG00000026959.14:B3", ("INH_medium_spiny", "INH_midbrain")),
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starsolo-root", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--primer-pairs", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--ggsashimi", type=Path)
    parser.add_argument("--run-plots", action="store_true")
    args = parser.parse_args()
    catalog = pd.read_csv(args.event_catalog, sep="\t", dtype=str)
    rows, events = {}, []
    for gene, event_id, cell_types in EVENTS:
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
    selected_cell_types = {cell_type for _, _, cell_types in EVENTS for cell_type in cell_types}
    barcode_groups, cell_sizes = read_primer_cell_groups(args.metadata, args.primer_pairs, selected_cell_types)
    bam_paths = sorted(args.starsolo_root.glob("*/spliced_pseudobulk.bam"))
    if not bam_paths:
        raise FileNotFoundError(f"no spliced pseudobulk BAMs under {args.starsolo_root}")
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(collect_bam_event_support, bam_paths, repeat(barcode_groups), repeat(tuple(events))))
    exon_blocks, junctions, event_umis = merge_support(results)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell_sizes.to_csv(args.output_dir / "cell_counts.tsv", sep="\t", index=False)
    support_rows = []
    for key, count in sorted(event_umis.items()):
        event_id, subject, cell_type, primer = key
        support_rows.append({"event": event_id, "subject": subject, "cell_type": cell_type, "primer": primer, "spliced_umis": count})
    pd.DataFrame(support_rows).to_csv(args.output_dir / "spliced_umi_support.tsv", sep="\t", index=False)
    manifest = {"normalization": "UMI-deduplicated spliced alignments per 1000 primer-specific half-cells", "runs": len(bam_paths), "events": []}
    palette = args.output_dir / "palette.tsv"
    palette.write_text("#1B9E77\n#D95F02\n#1B9E77\n#D95F02\n")
    event_ids = {gene: event_id for gene, event_id, _ in EVENTS}
    for event in events:
        event_dir = args.output_dir / event.event_id
        cell_types = rows[event.event_id + "_cell_types"]
        intron, exon, strand = write_ggsashimi_inputs(event_dir, event, exon_blocks, junctions, cell_sizes, cell_types)
        gtf = event_dir / f"{event.event_id}_paths.gtf"
        _event_gtf(gtf, event, rows[event.event_id])
        coordinates = f"{event.chromosome}:{event.start + 1}-{event.end + 1}"
        prefix = event_dir / f"{event.event_id}_priming_sashimi"
        manifest["events"].append({"gene": event.event_id, "test_id": event_ids[event.event_id], "cell_types": list(cell_types), "coordinates": coordinates, "intron_counts": str(intron.relative_to(args.output_dir)), "exon_counts": str(exon.relative_to(args.output_dir)), "path_annotation": str(gtf.relative_to(args.output_dir)), "figure": str(prefix.with_suffix('.pdf').relative_to(args.output_dir))})
        if args.run_plots:
            if args.ggsashimi is None:
                raise ValueError("--ggsashimi is required with --run-plots")
            subprocess.run(["python", str(args.ggsashimi), "--intron", str(intron), "--exon", str(exon), "--strand_info", str(strand), "--data_type", "sc", "--aggregation", "False", "-c", coordinates, "-g", str(gtf), "-o", str(prefix), "-M", "1", "--shrink", "--fix-y-scale", "--height", "0.8", "--ann-height", "0.7", "--width", "9", "--base-size", "11", "-P", str(palette), "-F", "pdf"], check=True)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
