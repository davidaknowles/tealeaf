#!/usr/bin/env python3
"""Describe frozen comparator rank-prefix representation, never select events."""

import argparse
import json
from pathlib import Path

import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash


def catalog_retention(ranked, reference, definitions):
    """Event presence is an upper bound, not proof a contrast was fitted."""
    if reference.feature_id.duplicated().any() or definitions.feature_id.duplicated().any():
        raise ValueError("unique frozen event families required")
    table = ranked.copy()
    old_events = set(reference.feature_id)
    old_genes = {str(value).split(".", 1)[0] for value in reference.gene_id}
    table["stable_gene"] = table.feature_id.str.removeprefix("SUPPA2:").str.split(";", n=1).str[0].str.split(".", n=1).str[0]
    table["prior_definition_supported"] = table.feature_id.isin(old_events)
    table["prior_gene_supported"] = table.stable_gene.isin(old_genes)
    table["full_catalog_geometry_status"] = table.feature_id.map(definitions.set_index("feature_id").status).fillna("absent from source catalog")
    table["representation"] = "gene absent from prior coverage/support family"
    table.loc[table.prior_gene_supported, "representation"] = "gene present, event absent from prior support family"
    table.loc[table.prior_definition_supported, "representation"] = "event present in prior support family"
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published-ranks", type=Path, required=True)
    parser.add_argument("--reference-family", type=Path, required=True)
    parser.add_argument("--full-definitions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new diagnostic output required")
    ranked = pd.read_csv(args.published_ranks, sep="\t", low_memory=False)
    ranked = ranked.loc[ranked.method.eq("SUPPA2 (full data)")].copy()
    if len(ranked) != 200 or ranked['rank'].duplicated().any() or set(ranked['rank']) != set(range(1, 201)):
        raise ValueError("unchanged complete published native top-200 prefix required")
    reference = pd.read_csv(args.reference_family, sep="\t")
    reference = reference.loc[reference.source.eq("original_binary")]
    definitions = pd.read_csv(args.full_definitions, sep="\t")
    table = catalog_retention(ranked, reference, definitions)
    records = []
    for cutoff in (100, 200):
        prefix = table.loc[table['rank'].le(cutoff)]
        for label, local in prefix.groupby("representation"):
            records.append(dict(cutoff=cutoff, representation=label, requested=len(local), unique_genes=local.stable_gene.nunique(), unique_event_definitions=local.feature_id.nunique(), marker_geometry_available=int(local.full_catalog_geometry_status.eq("ok").sum())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "rank_prefix.tsv.gz", sep="\t", index=False)
    summary = pd.DataFrame(records)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(input_hashes={str(path): file_hash(path) for path in (args.published_ranks, args.reference_family, args.full_definitions)}, scope="descriptive representation of unchanged published comparator prefixes only; new full catalog was frozen independently before this audit", limitations="Presence in the older supported event family does not prove that specific contrast was requested, fitted or LR-evaluable. Gene coverage, EC dimension and transcript support remain confounded in absence categories. No new significance ranking, cutoff, event filter, power or replication claim.", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
