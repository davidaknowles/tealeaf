"""Separate missing event expression from structural EC event nonidentifiability."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, event_path_index, supported_gene_transcripts
from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, local_test_design, modeled_gene_umis
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.event_paths import binary_event_information, mixture_event_information


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--null-panel", type=Path, required=True)
    parser.add_argument("--native-support", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with args.data_cache.open("rb") as handle:
        _, _, genes, gene_tx, gene_ecs, designs = pickle.load(handle)
    features = args.features.read_text().splitlines()
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    native = pd.read_csv(args.native_support, sep="\t")
    null_settings = json.loads(args.null_panel.read_text())
    metadata, counts, _, _, _, _ = filtered_inputs(args.data_cache, null_settings["candidate_settings"])
    screening = tuple(matrix.tocsc() for matrix in counts)
    contexts = {}
    records = [{"scope": "fixed random whole-family null panel", "test_id": key, "feature_id": key.rsplit("|", 3)[0]} for key in null_settings["requested_ids"]]
    records += [{"scope": "native top100 diagnostic, not a model-selection family", "test_id": row.feature_id + "|" + row.contrast_id, "feature_id": row.feature_id, "native_rank": row.native_rank} for row in native.loc[native.native_rank.le(100)].itertuples()]
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    cache, output = {}, []
    for record in records:
        feature = record["feature_id"]
        if feature not in cache:
            event = catalog.loc[feature]
            gene = lookup[canonical(event.gene_id)]
            transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
            paths = event_path_index(transcripts, features, event.included, event.excluded)
            values = {"gene_id": event.gene_id, "event_type": event.event_type, "n_transcripts": len(transcripts), "n_ecs": len(gene_ecs[gene]), "event_complete": paths is not None}
            if paths is not None:
                maps = tuple(np.asarray(mapping[gene_ecs[gene]][:, transcripts].toarray()) for mapping in designs)
                information = [binary_event_information(maps, paths, inclusion=psi) for psi in (.2, .5, .8)]
                values.update({"free_information_max": max(row["free_transcript_information"] for row in information), "fixed_information_max": max(row["fixed_mixture_information"] for row in information)})
                values["free_information_positive"] = values["free_information_max"] > 1e-6
                values["fixed_information_positive"] = values["fixed_information_max"] > 1e-6
                values["information_requires_fixed_shares"] = values["fixed_information_positive"] and not values["free_information_positive"]
                values["information_ratio"] = values["free_information_max"] / values["fixed_information_max"] if values["fixed_information_positive"] else np.nan
            cache[feature] = values
        detail = {}
        if record["scope"] == "fixed random whole-family null panel":
            event = catalog.loc[feature]
            gene = lookup[canonical(event.gene_id)]
            transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
            if gene not in contexts:
                settings = null_settings["candidate_settings"]
                umis = modeled_gene_umis(screening, designs, gene_ecs[gene], transcripts)
                specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=settings["min_gene_umis"], min_samples=settings["min_gene_samples"], min_celltype_mice=settings["min_celltype_mice"])
                contexts[gene] = {tuple(spec[0][-1]): spec[0][0] for spec in specs}
            levels = tuple(record["test_id"].rsplit("|", 2)[1:])
            rows = contexts[gene][levels]
            local, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(rows), 1)), local.mouse.to_numpy(), drop_zero=False)
            baseline = pooled_isoform_weights(base)
            paths = event_path_index(transcripts, features, event.included, event.excluded)
            if paths is not None:
                mass = baseline[paths >= 0].sum()
                totals = np.asarray([[matrix[labels == label].sum() for label in (0, 1)] for matrix in base.counts])
                fitted_info = mixture_event_information(base.compatibility, paths, np.maximum(baseline, 1e-12), totals)
                detail = {"pooled_event_mass": mass, "pooled_inclusion": baseline[paths == 0].sum() / mass, "pooled_gene_molecules": totals.sum(), "pooled_free_information": fitted_info["free_transcript_information"], "pooled_fixed_information": fitted_info["fixed_mixture_information"]}
        output.append({**record, **cache[feature], **detail})
    table = pd.DataFrame(output)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "events.tsv", sep="\t", index=False, na_rep="NA")
    summary = table.groupby("scope").agg(n_associations=("test_id", "size"), n_events=("feature_id", "nunique"), n_complete=("event_complete", "sum"), free_informative=("free_information_positive", "sum"), fixed_informative=("fixed_information_positive", "sum"), requires_fixed_shares=("information_requires_fixed_shares", "sum")).reset_index()
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({"source": str(args.data_cache), "interior_inclusion": [.2, .5, .8], "event_mass_if_outside": .7, "total_per_type_primer": 10000, "positive_information_threshold": 1e-6, "null": "two types with identical positive class-uniform transcript mixture, expected counts only", "profile": "type-specific within-class transcript shares and outside mass versus class-fixed shares", "limitation": "generic interior likelihood information, not a calibration result or actual expression-level eligibility; finite precision and annotation/mapping assumptions remain", "selection": "fixed null-family random panel plus labeled native-top100 diagnosis, not production/LR-selected model training", "production_changes": False}, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
