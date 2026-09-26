#!/usr/bin/env python3
"""Run SUPPA2 event quantification with a vectorized classical test.

SUPPA2 v2.4 generates the local-event catalogue and PSI values. Its classical
paired test is a Wilcoxon signed-rank test for each event followed by BH
correction. The upstream implementation loops over events in Python, so this
driver evaluates the same paired Wilcoxon normal approximation vectorized over
events after native `generateEvents` and `psiPerEvent`. The output keeps raw
p-values and reports the BH q-value separately.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import ndtr


EVENT_TYPES = ("SE", "A3", "A5", "MX", "RI")
EVENT_LABELS = {"A3": "A3SS", "A5": "A5SS", "MX": "MXE"}


def bh(values):
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if not finite.any():
        return out
    order = np.argsort(values[finite], kind="stable")
    ordered = values[finite][order]
    adjusted = np.minimum.accumulate((ordered * len(ordered) / np.arange(1, len(ordered) + 1))[::-1])[::-1]
    restored = np.empty_like(ordered)
    restored[order] = np.minimum(adjusted, 1.0)
    out[finite] = restored
    return out


def fast_paired_wilcoxon(differences, valid):
    """Vectorized normal-approximation paired Wilcoxon p-values.

    SUPPA2's classical test uses scipy's paired Wilcoxon test event by event.
    The rank statistic has a simple vectorized approximation, which avoids
    millions of Python-level scipy calls while retaining the same test family.
    """
    n_events = differences.shape[0]
    p_values = np.full(n_events, np.nan, dtype=float)
    nonzero = valid & (differences != 0)
    n = nonzero.sum(axis=1)
    selected = np.flatnonzero(n >= 1)
    if not len(selected):
        return p_values
    values = np.abs(differences[selected]).copy()
    mask = nonzero[selected]
    values[~mask] = np.inf
    order = np.argsort(values, axis=1, kind="stable")
    ranks = np.empty_like(values, dtype=float)
    positions = np.broadcast_to(np.arange(1, values.shape[1] + 1), values.shape)
    np.put_along_axis(ranks, order, positions, axis=1)
    ranks[~mask] = 0.0
    w_plus = np.sum(np.where((differences[selected] > 0) & mask, ranks, 0.0), axis=1)
    n_selected = n[selected]
    expected = n_selected * (n_selected + 1) / 4.0
    variance = n_selected * (n_selected + 1) * (2 * n_selected + 1) / 24.0
    z = (w_plus - expected) / np.sqrt(variance)
    p_values[selected] = 2.0 * ndtr(-np.abs(z))
    return p_values


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--columns", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--suppa-root", required=True, type=Path)
    parser.add_argument("--contrasts", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--event-output", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--minimum-pairs", type=int, default=8)
    return parser.parse_args()


def suppa_command(root, *args):
    return ["python3", str(root / "suppa.py"), *map(str, args)]


def load_events(work_dir):
    rows = []
    for event_type in EVENT_TYPES:
        path = work_dir / f"mouse_{event_type}_strict.ioe"
        table = pd.read_csv(path, sep="\t", dtype=str)
        for record in table.itertuples(index=False):
            event_id = str(record.event_id)
            included = [x for x in str(record.alternative_transcripts).split(",") if x]
            total = [x for x in str(record.total_transcripts).split(",") if x]
            excluded = sorted(set(total).difference(included))
            if not included or not excluded:
                continue
            rows.append({"feature_id": f"SUPPA2:{event_id}", "event_id": event_id, "event_type": EVENT_LABELS.get(event_type, event_type), "gene_id": str(record.gene_id), "gene_name": "", "included": ",".join(included), "excluded": ",".join(excluded)})
    return pd.DataFrame(rows).drop_duplicates("feature_id").reset_index(drop=True)


def write_expression(matrix, transcripts, row_lookup, samples, path):
    selected = [sample for sample in samples if sample in row_lookup]
    if not selected:
        return []
    values = matrix[[row_lookup[sample] for sample in selected]].toarray().T
    with path.open("w") as handle:
        handle.write("\t".join(selected) + "\n")
        for transcript, row in zip(transcripts, values):
            handle.write(transcript + "\t" + "\t".join(f"{value:.8g}" for value in row) + "\n")
    return selected


def run_native_psi(root, ioe, expression, output_prefix, log_path):
    psi_path = Path(str(output_prefix) + ".psi")
    if psi_path.exists():
        return psi_path
    command = suppa_command(root, "psiPerEvent", "-i", ioe, "-e", expression, "-o", output_prefix, "-m", "ERROR")
    with log_path.open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=log)
    return psi_path


def load_psi(path):
    table = pd.read_csv(path, sep="\t", index_col=0, na_values=["NA", "nan"])
    table.index = table.index.astype(str)
    return table


def main():
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    manifest_groups = [json.loads(path.read_text()) for path in args.contrasts]
    contrasts = [record for group in manifest_groups for record in group]
    fold_by_contrast = {record["contrast_id"]: fold for fold, group in enumerate(manifest_groups) for record in group}
    samples = sorted({sample for record in contrasts for sample in record["samples_a"] + record["samples_b"]})
    matrix = sparse.load_npz(args.matrix).tocsr()
    transcripts = [line.strip() for line in args.columns.read_text().splitlines()]
    rows = [line.strip() for line in args.rows.read_text().splitlines()]
    row_lookup = {row.split("__")[-1] + "__" + row.split("__")[0]: index for index, row in enumerate(rows)}
    cell_types = sorted({sample.rsplit("__", 1)[1] for sample in samples})
    ioe_prefix = args.work_dir / "mouse"
    ioe = args.work_dir / "mouse.ioe"
    if not ioe.exists():
        subprocess.run(suppa_command(args.suppa_root, "generateEvents", "-i", args.gtf, "-o", ioe_prefix, "-f", "ioe", "-e", "SE", "SS", "MX", "RI", "-m", "ERROR"), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        parts = [args.work_dir / f"mouse_{event_type}_strict.ioe" for event_type in EVENT_TYPES]
        with ioe.open("w") as handle:
            for index, part in enumerate(parts):
                lines = part.read_text().splitlines()
                if index:
                    lines = lines[1:]
                handle.write("\n".join(lines) + "\n")
    events = load_events(args.work_dir)
    args.event_output.parent.mkdir(parents=True, exist_ok=True)
    events.to_csv(args.event_output, sep="\t", index=False, compression="gzip")
    psi_by_cell = {}
    for cell_type in cell_types:
        expression = args.work_dir / f"{cell_type}.tsv"
        selected = write_expression(matrix, transcripts, row_lookup, [sample for sample in samples if sample.endswith("__" + cell_type)], expression)
        if not selected:
            continue
        psi_path = run_native_psi(args.suppa_root, ioe, expression, args.work_dir / f"{cell_type}.psi", args.work_dir / f"{cell_type}.log")
        psi_table = load_psi(psi_path)
        # SUPPA2 can omit an expression column that has no usable transcript
        # values, so use the columns actually returned by psiPerEvent.
        psi_by_cell[cell_type] = (list(psi_table.columns), psi_table)
    event_ids = events.event_id.to_numpy()
    result_rows = []
    for contrast in contrasts:
        level_a, level_b = contrast["level_a"], contrast["level_b"]
        if level_a not in psi_by_cell or level_b not in psi_by_cell:
            continue
        samples_a, psi_a = psi_by_cell[level_a]
        samples_b, psi_b = psi_by_cell[level_b]
        pairs = [(a, b) for a, b in zip(contrast["samples_a"], contrast["samples_b"]) if a in samples_a and b in samples_b]
        if len(pairs) < args.minimum_pairs:
            continue
        columns_a = [samples_a.index(a) for a, _ in pairs]
        columns_b = [samples_b.index(b) for _, b in pairs]
        first = psi_a.reindex(event_ids).to_numpy(dtype=float)[:, columns_a]
        second = psi_b.reindex(event_ids).to_numpy(dtype=float)[:, columns_b]
        valid = np.isfinite(first) & np.isfinite(second)
        enough = valid.sum(axis=1) >= args.minimum_pairs
        differences = second - first
        sums = np.nansum(np.where(valid, differences, np.nan), axis=1)
        counts = valid.sum(axis=1)
        effects = np.divide(
            sums,
            counts,
            out=np.full(len(events), np.nan, dtype=float),
            where=counts > 0,
        )
        p_values = np.full(len(events), np.nan)
        selected = np.flatnonzero(enough)
        if len(selected):
            p_values = fast_paired_wilcoxon(differences, valid)
        q_values = bh(p_values)
        for index in selected:
            event = events.iloc[index]
            result_rows.append({"method": "SUPPA2 (full data)", "contrast_id": contrast["contrast_id"], "effect": "cell_type", "stratum": contrast.get("stratum", "all"), "level_a": level_a, "level_b": level_b, "feature_id": event.feature_id, "event_type": event.event_type, "event_id": event.event_id, "p_value": float(p_values[index]), "q_value": float(q_values[index]), "effect_size": float(effects[index]), "n_subjects": int(valid[index].sum()), "fold": fold_by_contrast[contrast["contrast_id"]], "gene_id": event.gene_id, "gene_name": event.gene_name, "significant": bool(q_values[index] < 0.05), "criterion": "SUPPA2 classical paired Wilcoxon, BH q < 0.05"})
    result = pd.DataFrame(result_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(result):,} SUPPA2 full-data tests and {len(events):,} events")


if __name__ == "__main__":
    main()
