"""Retain all declared sequence-control cases and expose unsupported counts."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def summarize(table):
    ok = table.status.eq("ok")
    labels = table.full_positive_support.fillna(False).astype(str).str.lower()
    if not labels.isin(("true", "false")).all():
        raise ValueError("invalid full-support boolean")
    full = ok & labels.eq("true")
    paired = ok & table.sequence_KL_lower_bound.notna() & table.baseline_same_rows_KL_lower_bound.notna()
    result = dict(requested_gene_primer_rows=len(table), genes=table.gene_id.nunique(), successful_kernel_rows=int(ok.sum()), full_positive_support_rows=int(full.sum()), partial_positive_support_rows=int((ok & ~full & table.retained_molecule_fraction.gt(0)).sum()), median_retained_molecule_fraction=float(table.retained_molecule_fraction.median()), pooled_retained_molecule_fraction=float((table.input_molecules * table.retained_molecule_fraction.fillna(0)).sum() / table.input_molecules.sum()), paired_KL_rows=int(paired.sum()))
    for name, selected in (("supported_rows_only", paired), ("full_positive_support", paired & full)):
        result[f"{name}_median_sequence_KL_lower_bound"] = float(table.loc[selected, "sequence_KL_lower_bound"].median())
        result[f"{name}_median_baseline_KL_lower_bound"] = float(table.loc[selected, "baseline_same_rows_KL_lower_bound"].median())
        difference = table.loc[selected, "sequence_KL_lower_bound"] - table.loc[selected, "baseline_same_rows_KL_lower_bound"]
        result[f"{name}_median_KL_difference"] = float(difference.median())
        result[f"{name}_fraction_lower_KL"] = float(difference.lt(0).mean()) if len(difference) else np.nan
        result[f"{name}_median_sequence_achieved_KL"] = float(table.loc[selected, "sequence_KL"].median())
        result[f"{name}_median_baseline_achieved_KL"] = float(table.loc[selected, "baseline_same_rows_KL"].median())
        result[f"{name}_median_sequence_KL_certificate_gap"] = float((table.loc[selected, "sequence_KL"] - table.loc[selected, "sequence_KL_lower_bound"]).median())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=8)
    args = parser.parse_args()
    frames, metadata, summaries = [], [], []
    for source in ("parsimony_binary", "original_binary"):
        paths = [args.input_root / source / f"shard_{index}" for index in range(args.shard_count)]
        settings = [json.loads((path / "settings.json").read_text()) for path in paths]
        if any(value != settings[0] for value in settings) or settings[0]["shard_count"] != args.shard_count:
            raise ValueError("inconsistent sequence-control settings")
        pieces = [pd.read_csv(path / "genes.tsv.gz", sep="\t", dtype={"end_window": str}) for path in paths]
        table = pd.concat([piece for piece in pieces if len(piece)], ignore_index=True)
        observed = set(map(tuple, table[["gene_id", "read_length", "end_window", "primer"]].to_numpy()))
        expected = {(gene, length, str(window), primer) for gene in settings[0]["requested_genes"] for length in settings[0]["read_lengths"] for window in settings[0]["terminal_start_windows"] for primer in (0, 1)}
        if len(table) != len(observed) or observed != expected or not table.source.eq(source).all():
            raise ValueError("partial, duplicate or wrong-source sequence-control family")
        for field in ("retained_molecule_fraction", "exact_class_fraction"):
            if field not in table:
                table[field] = np.nan
            if not table[field].dropna().between(0, 1).all():
                raise ValueError("invalid sequence-control support fraction")
        for field in ("sequence_KL", "baseline_same_rows_KL", "sequence_KL_lower_bound", "baseline_same_rows_KL_lower_bound", "full_positive_support"):
            if field not in table:
                table[field] = np.nan
        if "background_fraction" in table and not table.background_fraction.eq(settings[0].get("background_fraction", 0.)).all():
            raise ValueError("mixed background recipes within sequence-control family")
        metadata.append(settings[0])
        frames.append(table)
        for (length, window, primer), group in table.groupby(["read_length", "end_window", "primer"]):
            summaries.append(dict(source=source, read_length=length, end_window=window, primer=primer, **summarize(group)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(args.output_dir / "genes.tsv.gz", sep="\t", index=False)
    result = pd.DataFrame(summaries)
    result.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(source_settings=metadata, inference="diagnostic only, no production maps or statistical procedure replaced", denominators="all declared genes, lengths, windows and primers including sequence/class failures; omitted molecule mass remains explicit", KL="sequence and baseline maps evaluated on exactly the same represented EC rows; partial-row discrepancy is not full-model validation", caution="gene-local error-free single-read matching need not reproduce the real UMI EC construction, primer distribution or transcriptome-wide ambiguity; no split/own-ranked LR claim"), indent=2) + "\n")
    print(result.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
