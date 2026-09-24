#!/usr/bin/env python3
"""Run a transcript-PSI, SUPPA-style classical event comparison.

SUPPA itself is not installed on the analysis system.  This implementation
uses the same transcript-usage input as Tealeaf, constructs the standard
annotated SE/A3SS/A5SS/MXE/RI event classes, and applies SUPPA's classical
paired PSI test (paired t test followed by within-contrast BH correction). It
is therefore labelled ``SUPPA transcript PSI`` rather than implying that a
different quantifier was used.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_rel


ATTR = re.compile(r'(transcript_id|gene_id|gene_name) "([^"]+)"')
EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")


def bh(values):
    values = np.asarray(values, dtype=float)
    output = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if not finite.any():
        return output
    order = np.argsort(values[finite])
    ordered = values[finite][order]
    adjusted = np.minimum.accumulate((ordered * len(ordered) / np.arange(1, len(ordered) + 1))[::-1])[::-1]
    result = np.empty_like(ordered)
    result[order] = np.minimum(adjusted, 1.0)
    output[finite] = result
    return output


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--columns", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--events", required=True, type=Path, help="rMATS fromGTF directory used as the standard event catalogue")
    parser.add_argument("--contrasts", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--event-output", required=True, type=Path)
    parser.add_argument("--minimum-pairs", type=int, default=8)
    return parser.parse_args()


def interval(start, end):
    return int(start), int(end)


def load_transcripts(gtf):
    exons = defaultdict(list)
    gene = {}
    name = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(ATTR.findall(fields[8]))
            transcript = attrs.get("transcript_id")
            if not transcript:
                continue
            transcript = transcript.split(".")[0]
            exons[transcript].append((fields[0], fields[6], int(fields[3]) - 1, int(fields[4])))
            gene[transcript] = attrs.get("gene_id", "").split(".")[0]
            name[transcript] = attrs.get("gene_name", "")
    normalized = {}
    for transcript, records in exons.items():
        normalized[transcript] = sorted(set(records), key=lambda x: x[2])
    return normalized, gene, name


def build_indexes(transcript_exons):
    exact = defaultdict(set)
    adjacent = defaultdict(set)
    for transcript, exons in transcript_exons.items():
        local = [(c, st, s, e) for c, st, s, e in exons]
        for chromosome, strand, start, end in local:
            exact[(chromosome, strand, start, end)].add(transcript)
        for first, second in zip(local, local[1:]):
            if first[0] == second[0] and first[1] == second[1]:
                adjacent[(first[0], first[1], first[3], second[2])].add(transcript)
    return exact, adjacent


def add_event(event_rows, event_type, row, transcript_exons, transcript_gene, exact_index, adjacent):
    chromosome = str(row["chr"]).strip('"')
    strand = str(row["strand"]).strip('"')
    def coord(column_start, column_end):
        return chromosome, strand, int(row[column_start]), int(row[column_end])

    if event_type == "SE":
        middle = coord("exonStart_0base", "exonEnd")
        upstream = coord("upstreamES", "upstreamEE")
        downstream = coord("downstreamES", "downstreamEE")
        included = exact_index[middle] & exact_index[upstream] & exact_index[downstream]
        excluded = exact_index[upstream] & exact_index[downstream] & adjacent[(chromosome, strand, upstream[3], downstream[2])]
    elif event_type in ("A3SS", "A5SS"):
        long = coord("longExonStart_0base", "longExonEnd")
        short = coord("shortES", "shortEE")
        flank = coord("flankingES", "flankingEE")
        included = exact_index[long] & exact_index[flank]
        excluded = exact_index[short] & exact_index[flank]
    elif event_type == "MXE":
        first = coord("1stExonStart_0base", "1stExonEnd")
        second = coord("2ndExonStart_0base", "2ndExonEnd")
        upstream = coord("upstreamES", "upstreamEE")
        downstream = coord("downstreamES", "downstreamEE")
        included = exact_index[upstream] & exact_index[first] & exact_index[downstream]
        excluded = exact_index[upstream] & exact_index[second] & exact_index[downstream]
    else:
        retained = coord("riExonStart_0base", "riExonEnd")
        upstream = coord("upstreamES", "upstreamEE")
        downstream = coord("downstreamES", "downstreamEE")
        included = exact_index[retained]
        excluded = exact_index[upstream] & exact_index[downstream] & adjacent[(chromosome, strand, upstream[3], downstream[2])]
    included -= excluded
    included = sorted(included)
    excluded = sorted(excluded)
    if not included or not excluded:
        return None
    gene_ids = {transcript_gene.get(t, "") for t in included + excluded if transcript_gene.get(t)}
    return {"event_type": event_type, "event_id": str(row["ID"]), "gene_id": sorted(gene_ids)[0] if len(gene_ids) == 1 else "", "included": included, "excluded": excluded}


def event_catalog(event_dir, gtf):
    transcript_exons, transcript_gene, transcript_name = load_transcripts(gtf)
    exact, adjacent = build_indexes(transcript_exons)
    rows = []
    for event_type in EVENT_TYPES:
        path = event_dir / f"fromGTF.{event_type}.txt"
        table = pd.read_csv(path, sep="\t", low_memory=False)
        for row in table.to_dict("records"):
            event = add_event(rows, event_type, row, transcript_exons, transcript_gene, exact, adjacent)
            if event:
                event["gene_name"] = next((transcript_name.get(t, "") for t in event["included"] + event["excluded"] if transcript_name.get(t)), "")
                rows.append(event)
    return rows


def main():
    args = parse_args()
    matrix = sparse.load_npz(args.matrix).tocsr()
    row_names = args.rows.read_text().splitlines()
    transcript_names = args.columns.read_text().splitlines()
    # GTF and transcript matrices may differ only by version suffix.
    matrix_index = defaultdict(list)
    for index, transcript in enumerate(transcript_names):
        matrix_index[transcript.split(".")[0]].append(index)
    row_lookup = {}
    for index, label in enumerate(row_names):
        pieces = label.split("__")
        if len(pieces) >= 3:
            cell_type, condition, subject = "__".join(pieces[:-2]), pieces[-2], pieces[-1]
            row_lookup[f"{subject}__{cell_type}"] = index
    events = event_catalog(args.events, args.gtf)
    args.event_output.parent.mkdir(parents=True, exist_ok=True)
    event_table = pd.DataFrame([{**event, "included": ",".join(event["included"]), "excluded": ",".join(event["excluded"]), "feature_id": f"SUPPA:{event['event_type']}:{event['event_id']}"} for event in events])
    event_table.to_csv(args.event_output, sep="\t", index=False, compression="gzip")
    inc_rows, inc_cols, exc_rows, exc_cols = [], [], [], []
    for event_index, event in enumerate(events):
        for transcript in event["included"]:
            inc_rows.extend([event_index] * len(matrix_index.get(transcript, [])))
            inc_cols.extend(matrix_index.get(transcript, []))
        for transcript in event["excluded"]:
            exc_rows.extend([event_index] * len(matrix_index.get(transcript, [])))
            exc_cols.extend(matrix_index.get(transcript, []))
    shape = (len(events), len(transcript_names))
    inc_map = sparse.csr_matrix((np.ones(len(inc_rows)), (inc_rows, inc_cols)), shape=shape)
    exc_map = sparse.csr_matrix((np.ones(len(exc_rows)), (exc_rows, exc_cols)), shape=shape)
    transcript_values = matrix.T.tocsr()
    inc_values = inc_map @ transcript_values
    exc_values = exc_map @ transcript_values
    psi = np.divide(inc_values.toarray(), (inc_values + exc_values).toarray(), out=np.full((len(events), matrix.shape[0]), np.nan), where=(inc_values + exc_values).toarray() > 0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_handle = gzip.open(args.output, "wt")
    wrote_header = False
    for fold, path in enumerate(args.contrasts):
        contrasts = json.loads(path.read_text())
        for contrast in contrasts:
            a = [row_lookup.get(sample) for sample in contrast["samples_a"]]
            b = [row_lookup.get(sample) for sample in contrast["samples_b"]]
            pairs = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
            if len(pairs) < args.minimum_pairs:
                continue
            ai = np.asarray([x for x, _ in pairs])
            bi = np.asarray([y for _, y in pairs])
            first = psi[:, ai]
            second = psi[:, bi]
            valid = np.isfinite(first) & np.isfinite(second)
            enough = valid.sum(axis=1) >= args.minimum_pairs
            deltas = np.where(valid, second - first, np.nan)
            p_values = np.full(len(events), np.nan)
            effects = np.divide(np.nansum(deltas, axis=1), valid.sum(axis=1), out=np.full(len(events), np.nan), where=valid.sum(axis=1) > 0)
            selected_indices = np.flatnonzero(enough)
            if len(selected_indices):
                p_values[selected_indices] = ttest_rel(second[selected_indices], first[selected_indices], axis=1, nan_policy="omit").pvalue
            selected_rows = []
            for event_index in selected_indices:
                event = events[event_index]
                selected_rows.append({"method": "SUPPA transcript PSI", "contrast_id": contrast["contrast_id"], "effect": "cell_type", "stratum": contrast.get("stratum", "all"), "level_a": contrast["level_a"], "level_b": contrast["level_b"], "feature_id": f"SUPPA:{event['event_type']}:{event['event_id']}", "event_type": event["event_type"], "event_id": event["event_id"], "p_value": float(p_values[event_index]), "effect_size": float(effects[event_index]), "n_subjects": int(valid[event_index].sum()), "fold": fold, "gene_id": event["gene_id"], "gene_name": event["gene_name"]})
            result = pd.DataFrame(selected_rows)
            if len(result):
                result["q_value"] = bh(result["p_value"])
                result["significant"] = result["q_value"] < 0.05
                result["criterion"] = "paired transcript PSI t test, BH q < 0.05"
                result.to_csv(output_handle, sep="\t", index=False, header=not wrote_header)
                wrote_header = True
    output_handle.close()
    print(f"wrote SUPPA transcript-PSI tests and {len(events):,} events")


if __name__ == "__main__":
    main()
