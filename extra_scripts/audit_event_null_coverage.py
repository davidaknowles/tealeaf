"""Frozen-coverage count-null diagnostics including every failed hypothesis."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, local_test_design, modeled_gene_umis
from extra_scripts.summarize_ec_count_null import validate_requested_trials
from tealeaf.sc.replication_audit import coverage_correlation


def attach_frozen_coverage(trials, attributes):
    """Attach one pre-fit attribute vector to each trial, without case loss."""
    required = ("test_id", "median_gene_umis", "expected_subjects", "transcripts", "ecs")
    if any(column not in attributes for column in required) or attributes.test_id.duplicated().any() or set(trials.test_id) != set(attributes.test_id):
        raise ValueError("frozen coverage must include every unique requested hypothesis")
    if not np.isfinite(attributes[list(required[1:])].to_numpy(float)).all() or not attributes.median_gene_umis.gt(0).all() or not attributes.expected_subjects.ge(4).all():
        raise ValueError("positive pre-fit coverage and at least four eligible subjects required")
    result = trials.merge(attributes, on="test_id", how="left", validate="many_to_one", sort=False)
    if len(result) != len(trials):
        raise ValueError("coverage attachment changed the requested trial family")
    labels = result.converged.astype(str).str.lower()
    if not labels.isin(("true", "false")).all():
        raise ValueError("invalid null-fit availability")
    result["converged"] = labels.eq("true")
    if not result.loc[~result.converged, ["p_value", "raw_p_value"]].eq(1).all().all():
        raise ValueError("failed null trials must remain at p1")
    for column in ("p_value", "raw_p_value"):
        if not np.isfinite(result[column]).all() or not result[column].between(0, 1).all():
            raise ValueError("invalid count-null p-value")
    return result


def coverage_summaries(table):
    """All-trial rejection denominators, and separately labeled fit-only rho."""
    strata, correlations = [], []
    for (scenario, strategy, quartile), local in table.groupby(["scenario", "strategy", "coverage_quartile"], observed=True):
        row = dict(scenario=scenario, strategy=strategy, coverage_quartile=quartile, requested_trials=len(local), hypotheses=local.test_id.nunique(), fitted_trials=int(local.converged.sum()), median_gene_umis=float(local.median_gene_umis.median()), minimum_gene_umis=float(local.median_gene_umis.min()), maximum_gene_umis=float(local.median_gene_umis.max()))
        for name, column in (("native", "raw_p_value"), ("calibrated", "p_value")):
            for threshold in (.05, .01):
                count = int(local[column].le(threshold).sum())
                row[f"{name}_rejections_{threshold}"] = count
                row[f"{name}_rate_{threshold}"] = count / len(local)
        strata.append(row)
    for (scenario, strategy), local in table.groupby(["scenario", "strategy"], observed=True):
        for draw, selected in [("all repeated draws", local), *[(str(draw), group) for draw, group in local.groupby("draw")]]:
            for scope, frame in (("all requested trials", selected), ("successful-fit rho only, not calibration denominator", selected.loc[selected.converged])):
                controls = frame[["expected_subjects", "transcripts", "ecs"]]
                for name, column in (("native", "raw_p_value"), ("calibrated", "p_value")):
                    correlations.append(dict(scenario=scenario, strategy=strategy, draw=draw, scope=scope, p_scale=name, **coverage_correlation(frame[column], frame.median_gene_umis, controls)))
    return pd.DataFrame(strata), pd.DataFrame(correlations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True, help="Complete assessed strict/biological/nuisance null panels.")
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    panels, settings = [], []
    for scenario in ("strict", "biological", "nuisance"):
        root = args.input_root / scenario
        recipe = json.loads((root / "manifest.json").read_text())["settings"]
        panel = pd.read_csv(root / "tests.tsv.gz", sep="\t")
        validate_requested_trials(panel, recipe)
        panels.append(panel.assign(scenario=scenario))
        settings.append(recipe)
    identities = [set(value["requested_ids"]) for value in settings]
    if any(value != identities[0] for value in identities) or any(value["candidate_settings"] != settings[0]["candidate_settings"] for value in settings) or any(value["draws"] != settings[0]["draws"] for value in settings):
        raise ValueError("coverage comparison needs the same frozen hypotheses and screening across scenarios")
    expected_hash = settings[0].get("input_manifest_sha256")
    if expected_hash is not None:
        actual_hash = hashlib.sha256((args.data_cache.parent / "manifest.json").read_bytes()).hexdigest()
        if actual_hash != expected_hash or any(value.get("input_manifest_sha256") != expected_hash for value in settings):
            raise ValueError("coverage data cache differs from the declared count-null input")
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(args.data_cache, settings[0]["candidate_settings"])
    counts = tuple(value.tocsc() for value in counts)
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    candidate = settings[0]["candidate_settings"]
    attributes, contexts = [], {}
    for test_id in sorted(identities[0]):
        feature, factor, first, second = test_id.rsplit("|", 3)
        if factor != "cell_type" or not feature.startswith("SUPPA2:"):
            raise ValueError("unexpected frozen count-null hypothesis")
        gene_id = canonical(feature.removeprefix("SUPPA2:").split(";", 1)[0])
        gene = lookup[gene_id]
        if gene not in contexts:
            transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
            umis = modeled_gene_umis(counts, designs, gene_ecs[gene], transcripts)
            specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=candidate["min_gene_umis"], min_samples=candidate["min_gene_samples"], min_celltype_mice=candidate["min_celltype_mice"])
            contexts[gene] = umis, transcripts, {tuple(spec[0][-1]): spec[0][0] for spec in specs}
        umis, transcripts, rows_by_levels = contexts[gene]
        rows = rows_by_levels[(first, second)]
        local, _, labels = local_test_design(metadata, rows, (first, second), "cell_type_pairwise")
        subjects = local.mouse.astype(str).to_numpy()
        paired = set(subjects[labels == 0]) & set(subjects[labels == 1])
        paired_rows = np.asarray(rows)[np.isin(subjects, list(paired))]
        if len(paired_rows) != 2 * len(paired):
            raise ValueError("frozen paired coverage requires one row per subject/type")
        attributes.append(dict(test_id=test_id, gene_id=gene_id, median_gene_umis=float(np.median(umis[paired_rows])), total_gene_umis=float(umis[paired_rows].sum()), expected_subjects=len(paired), transcripts=len(transcripts), ecs=len(gene_ecs[gene])))
    attributes = pd.DataFrame(attributes)
    attributes["coverage_quartile"] = pd.qcut(attributes.median_gene_umis, 4, labels=False, duplicates="drop") + 1
    if attributes.coverage_quartile.isna().any():
        raise ValueError("frozen coverage has no usable quartile separation")
    table = attach_frozen_coverage(pd.concat(panels, ignore_index=True), attributes)
    strata, correlations = coverage_summaries(table)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    attributes.to_csv(args.output_dir / "frozen_hypothesis_coverage.tsv", sep="\t", index=False, na_rep="NA")
    strata.to_csv(args.output_dir / "coverage_quartiles.tsv", sep="\t", index=False, na_rep="NA")
    correlations.to_csv(args.output_dir / "coverage_correlations.tsv", sep="\t", index=False, na_rep="NA")
    hypotheses = table.groupby(["scenario", "strategy", "test_id"], observed=True).agg(requested_trials=("p_value", "size"), fitted_trials=("converged", "sum"), native_rejections_0_05=("raw_p_value", lambda values: int(values.le(.05).sum())), calibrated_rejections_0_05=("p_value", lambda values: int(values.le(.05).sum())), calibrated_rejections_0_01=("p_value", lambda values: int(values.le(.01).sum()))).reset_index().merge(attributes, on="test_id", how="left", validate="many_to_one")
    hypotheses.to_csv(args.output_dir / "hypothesis_coverage_summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(hypotheses=len(attributes), trials=len(table), selection="same frozen random hypotheses, including failures, no coverage or p-value selected panel", coverage="median original-input gene UMI total over the eligible paired subject/type rows, reconstructed independently of simulation fitting success", quartiles="defined on unique requested hypotheses before inspecting null p-values, shared across scenarios and draws", partial_rho_controls=["expected eligible subjects", "supported transcripts", "gene ECs"], dependence="four draws per hypothesis are not independent biological observations; rho is also reported per draw, no independent-binomial confidence claim", calibration_denominator="all requested trials, failed fits at p1; successful-fit-only rho is separately labeled and does not alter rejection denominators", data_cache=str(args.data_cache.resolve()), panel_settings=settings, production_changes=False), indent=2) + "\n")
    print(strata.to_string(index=False), flush=True)
    print(correlations.loc[correlations.draw.eq("all repeated draws")].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
