"""Observed-count score diagnostics, not a selected-family benchmark."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, local_test_design, modeled_gene_umis
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts, event_path_index
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.path_score_mixed import shared_path_score_components, aggregate_path_scores, paired_score_reporting, MODEL_VERSION
from tealeaf.sc.replication_audit import complete_cluster_fit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary"), required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--native-support", type=Path, required=True)
    parser.add_argument("--null-panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=8)
    parser.add_argument("--score-coordinate", choices=("ilr", "proportion"), default="ilr")
    parser.add_argument("--information-metric", choices=("absolute", "reference"), default="absolute")
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid diagnostic shard")
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    if settings["subject_fold"] is not None or settings["min_gene_umis"] != 25:
        raise ValueError("observed diagnostic requires the full-data production coverage recipe")
    source = args.cache / f"{args.source}_paired"
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(source / "prepared.pkl", settings)
    if metadata.duplicated(["mouse", "cell_type"]).any():
        raise ValueError("one pseudobulk per subject/type required")
    features = (source / "features.txt").read_text().splitlines()
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    native = pd.read_csv(args.native_support, sep="\t")
    records = [{"scope": "fixed random real-data diagnostic", "test_id": key, "feature_id": key.rsplit("|", 3)[0], "levels": tuple(key.rsplit("|", 2)[1:])} for key in json.loads(args.null_panel.read_text())["requested_ids"]]
    records += [{"scope": "native top100 diagnostic, not an inference family", "test_id": row.feature_id + "|" + row.contrast_id, "feature_id": row.feature_id, "levels": tuple(row.contrast_id.removeprefix("cell_type__").split("__")), "native_rank": row.native_rank, "long_read_effect": row.long_read_effect} for row in native.loc[native.native_rank.le(100)].itertuples()]
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    screening = tuple(value.tocsc() for value in counts)
    contexts, output, failures = {}, [], []
    for record in records[args.shard_index::args.shard_count]:
        result_row = {**record, "source": args.source, "p_value": 1., "converged": False}
        try:
            event = catalog.loc[record["feature_id"]]
            gene = lookup[canonical(event.gene_id)]
            if gene not in contexts:
                if len(gene_ecs[gene]) > settings["max_ecs"]:
                    raise ValueError("gene exceeds the frozen EC coverage ceiling")
                transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
                umis = modeled_gene_umis(screening, designs, gene_ecs[gene], transcripts)
                specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=settings["min_gene_umis"], min_samples=settings["min_gene_samples"], min_celltype_mice=settings["min_celltype_mice"])
                contexts[gene] = transcripts, {tuple(spec[0][-1]): spec[0][0] for spec in specs}
            transcripts, row_lookup = contexts[gene]
            rows = row_lookup[tuple(sorted(record["levels"]))]
            local, _, labels = local_test_design(metadata, rows, record["levels"], "cell_type_pairwise")
            subjects = local.mouse.to_numpy()
            paths = event_path_index(transcripts, features, event.included, event.excluded)
            if paths is None:
                raise ValueError("event lacks complete source transcript support")
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(rows), 1)), subjects, drop_zero=False)
            baseline = pooled_isoform_weights(base)
            components = shared_path_score_components(base, paths, labels, subjects, baseline=baseline, max_iter=2000, null_multistart=True, reporting_concentration=1., score_coordinate=args.score_coordinate)
            report = paired_score_reporting(components)
            result_row.update(n_expected_subjects=len(np.unique(subjects)), n_fitted_subjects=len(components.subject_ids), n_reported_subjects=report["n_reported_subjects"], report_complete=report["complete"], report_effect=report["effect"][0], pooled_inclusion=baseline[paths == 0].sum() / baseline[paths >= 0].sum(), pooled_event_mass=baseline[paths >= 0].sum(), n_interior_starts_selected=sum(fit.selected_start == 1 for fit in components.null_fits))
            result = aggregate_path_scores(components, information_metric=args.information_metric)
            complete = complete_cluster_fit(result, len(np.unique(subjects)))
            result_row.update(p_value=result["p_value"] if complete else 1., converged=complete, score_effect=result["mean_difference"][0], n_informative_subjects=result["n_subjects"], chi_square_p_value=result["chi_square_p_value"], statistic=result["statistic"])
            if "long_read_effect" in record:
                result_row.update(report_agrees=bool(report["effect"][0] * record["long_read_effect"] > 0) if report["complete"] else np.nan, score_agrees=bool(result["mean_difference"][0] * record["long_read_effect"] > 0) if complete else np.nan)
        except (ValueError, KeyError, np.linalg.LinAlgError) as error:
            result_row["error"] = repr(error)
            failures.append({"test_id": record["test_id"], "scope": record["scope"], "error": repr(error)})
        output.append(result_row)
        print(f"{record['test_id']}, complete={result_row['converged']}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output_dir / "tests.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "settings.json").write_text(json.dumps({"source": args.source, "candidate_settings": settings, "model_version": MODEL_VERSION, "score_coordinate": args.score_coordinate, "information_metric": args.information_metric, "max_iter": 2000, "null_multistart": True, "reporting_concentration": 1, "n_requested": len(records), "shard_count": args.shard_count, "selection": "fixed random null-family identifiers and labeled native top100 diagnosis; no selected-family BH or full-universe LR claim", "production_changes": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
