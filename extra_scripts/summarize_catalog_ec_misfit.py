"""Collate whole-family categorical-mixture diagnostics, not splicing p-values."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def summary(table):
    positive = table.positive_counts
    bounded = positive & table.integral_counts & table.log_probability_bound.notna()
    result = dict(gene_primer_rows=len(table), genes=table.gene_id.nunique(), positive_rows=int(positive.sum()), integral_positive_rows=int((positive & table.integral_counts).sum()), converged_positive_rows=int((positive & table.fit_converged).sum()), sampling_bound_rows=int(bounded.sum()))
    for field in ("molecules", "KL", "KL_lower_bound", "KL_dual_gap"):
        result[f"median_{field}"] = float(table.loc[positive, field].median())
    if "supported_molecule_fraction" in table:
        result["median_supported_molecule_fraction"] = float(table.supported_molecule_fraction.median())
        result["rows_with_mapping_support_loss"] = int(table.supported_molecule_fraction.lt(1 - 1e-12).sum())
    for cutoff in (.01, .001, 1e-6):
        rejected = bounded & table.log_probability_bound.lt(np.log(cutoff))
        result[f"incompatible_{cutoff}_rows"] = int(rejected.sum())
        result[f"incompatible_{cutoff}_fraction_all_rows"] = float(rejected.mean())
        result[f"incompatible_{cutoff}_fraction_bound_rows"] = float(rejected.sum() / bounded.sum()) if bounded.any() else np.nan
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sources", nargs="+", choices=("parsimony_binary", "original_binary", "prepared_weighted"), default=("parsimony_binary", "original_binary"))
    args = parser.parse_args()
    frames, summaries, metadata = [], [], []
    if len(set(args.sources)) != len(args.sources):
        raise ValueError("duplicate diagnostic sources")
    for source in args.sources:
        directories = [args.input_root / source / f"shard_{index}" for index in range(8)]
        settings = [json.loads((path / "settings.json").read_text()) for path in directories]
        if any(value != settings[0] for value in settings) or settings[0]["shard_count"] != 8:
            raise ValueError("inconsistent EC-model diagnostic settings")
        pieces = [pd.read_csv(path / "genes.tsv.gz", sep="\t") for path in directories]
        if not any(len(piece) for piece in pieces):
            raise ValueError("empty requested EC-model diagnostic family")
        table = pd.concat([piece for piece in pieces if len(piece)], ignore_index=True)
        observed = set(map(tuple, table[["gene_id", "primer"]].to_numpy()))
        expected = {(gene, primer) for gene in settings[0]["requested_genes"] for primer in (0, 1)}
        if len(table) != len(observed) or observed != expected or not table.source.eq(source).all():
            raise ValueError("incomplete, duplicate or wrong-source gene/primer family")
        for field in ("positive_counts", "integral_counts", "fit_converged"):
            normalized = table[field].astype(str).str.lower()
            if not normalized.isin(("true", "false")).all():
                raise ValueError("invalid diagnostic boolean field")
            table[field] = normalized.eq("true")
        if (table.KL_lower_bound.dropna() < 0).any() or (table.KL_lower_bound.dropna() > table.KL.dropna() + 1e-12).any() or table.loc[~table.integral_counts, "log_probability_bound"].notna().any():
            raise ValueError("invalid KL certificate or fractional-count probability bound")
        frames.append(table)
        metadata.append(settings[0])
        summaries.append(dict(source=source, primer="all", **summary(table)))
        for primer, group in table.groupby("primer"):
            summaries.append(dict(source=source, primer=str(primer), **summary(group)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(args.output_dir / "genes.tsv.gz", sep="\t", index=False)
    result = pd.DataFrame(summaries)
    result.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(source_settings=metadata, scope="complete screened gene family, per-primer unrestricted transcript mixtures", bound="conservative finite categorical model-compatibility reference bound using a certified KL lower bound, not a differential-splicing p-value", assumptions="fixed maps/categories and integral independent categorical molecules conditional on their observation-specific transcript mixtures; arbitrary mixtures and unequal observation totals allowed", limitations="maps may be estimated and EC retention is data-dependent, so this fixed-map reference bound is not a formal unconditional test on the selected real-data inputs; does not identify a unique cause, and molecule dependence, missing transcripts, read opportunities or primer/mapping bias can violate the working model; not proof of a different method's calibration or power", interpretation="diagnostic only, no gene exclusions, no revised discovery counts, no own-ranked LR or split superiority claim", production_changes=False), indent=2) + "\n")
    print(result.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
