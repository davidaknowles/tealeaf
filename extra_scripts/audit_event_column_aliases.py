"""Audit exact null-fitting redundancy on the whole screened event catalog."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts, event_path_index
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, modeled_gene_umis
from tealeaf.sc.path_aliases import exact_path_column_groups


def audit(cache, source, settings, catalog):
    control = cache / f"{source}_paired"
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(control / "prepared.pkl", settings)
    features = (control / "features.txt").read_text().splitlines()
    screening = tuple(value.tocsc() for value in counts)
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    rows = []
    for gene_id, events in catalog.groupby("canonical_gene", sort=False):
        gene = lookup.get(gene_id)
        if gene is None or not 0 < len(gene_ecs[gene]) <= settings["max_ecs"]:
            continue
        transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
        if len(transcripts) < 2:
            continue
        umis = modeled_gene_umis(screening, designs, gene_ecs[gene], transcripts)
        specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=settings["min_gene_umis"], min_samples=settings["min_gene_samples"], min_celltype_mice=settings["min_celltype_mice"])
        if not specs:
            continue
        mappings = tuple(value[gene_ecs[gene]][:, transcripts] for value in designs)
        for event in events.itertuples():
            paths = event_path_index(transcripts, features, event.included, event.excluded)
            if paths is None:
                continue
            groups = exact_path_column_groups(mappings, paths)
            rows.append(dict(source=source, gene_id=gene_id, feature_id=event.feature_id, transcripts=len(transcripts), unique_within_path_columns=len(groups), max_alias_group=max(map(len, groups)), n_ecs=len(gene_ecs[gene]), n_contrasts=len(specs), full_null_dimension=2 * len(transcripts) - 3, alias_null_dimension=2 * len(groups) - 3))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    if settings["subject_fold"] is not None:
        raise ValueError("use the frozen full-data coverage recipe")
    catalog = pd.read_csv(args.event_catalog, sep="\t")
    catalog["canonical_gene"] = catalog.gene_id.map(canonical)
    frames, summary = [], []
    for source in ("parsimony_binary", "original_binary"):
        table = audit(args.cache, source, settings, catalog)
        if table.empty:
            raise ValueError("whole catalog audit yielded no supported contexts")
        frames.append(table)
        weighted = table.n_contrasts.to_numpy(float)
        row = dict(source=source, genes=table.gene_id.nunique(), events=len(table), event_contrasts=int(weighted.sum()), events_with_alias=int(table.max_alias_group.gt(1).sum()), aliased_contrast_fraction=float(np.average(table.max_alias_group.gt(1), weights=weighted)), weighted_dimension_ratio=float(np.sum(weighted * table.alias_null_dimension) / np.sum(weighted * table.full_null_dimension)), maximum_transcripts=int(table.transcripts.max()), median_transcripts=float(table.transcripts.median()))
        summary.append(row)
        print(row, flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(args.output_dir / "events.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(candidate_settings=settings, selection="whole supported catalog and frozen full-data gene coverage, never significance or LR selection", equivalence="exact same column values in all primers, within the same event class", limitation="dimension diagnostic only, does not silently change nuisance prior multiplicities or optimizer constraints", production_changes=False), indent=2) + "\n")


if __name__ == "__main__":
    main()
