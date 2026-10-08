"""Whole supported-catalog path information, never a selected-hit benchmark."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts, event_path_index
from tealeaf.sc.ec_glmm import subset_gene_data
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.event_paths import binary_event_mixture, mixture_event_information


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary"), required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--alias-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=8)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid information-audit shard")
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    if settings["subject_fold"] is not None:
        raise ValueError("whole information audit requires the frozen full-data recipe")
    source = args.cache / f"{args.source}_paired"
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(source / "prepared.pkl", settings)
    counts = tuple(value.tocsc() for value in counts)
    features = (source / "features.txt").read_text().splitlines()
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    family = pd.read_csv(args.alias_audit, sep="\t")
    family = family.loc[family.source.eq(args.source)].copy()
    if family.empty or family.feature_id.duplicated().any():
        raise ValueError("empty or duplicate whole screened event family")
    gene_ids = sorted(family.gene_id.unique())
    selected = set(gene_ids[args.shard_index::args.shard_count])
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    output = []
    for gene_number, (gene_id, events) in enumerate(family.loc[family.gene_id.isin(selected)].groupby("gene_id", sort=False)):
        gene = lookup[canonical(gene_id)]
        transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
        base, _, _ = subset_gene_data(counts, designs, transcripts, gene_ecs[gene], np.ones((len(metadata), 1)), np.zeros(len(metadata)), drop_zero=False)
        baseline, converged = pooled_isoform_weights(base, return_status=True)
        # Preserve every event even if its pooled initialization is incomplete.
        baseline = np.maximum(baseline, 1e-12)
        baseline /= baseline.sum()
        primer_totals = np.asarray([value.sum() for value in base.counts])
        fractions = primer_totals / primer_totals.sum()
        interior_totals = np.tile(fractions[:, None], (1, 2))
        pooled_totals = np.tile(primer_totals[:, None] / 2, (1, 2))
        for event_record in events.itertuples():
            event = catalog.loc[event_record.feature_id]
            paths = event_path_index(transcripts, features, event.included, event.excluded)
            if paths is None or len(transcripts) != event_record.transcripts:
                raise ValueError("prepared input no longer matches the frozen audited family")
            weights = (paths == 0, paths == 1, paths < 0)
            assigned = np.zeros(2)
            for mapping, total in zip(base.compatibility, primer_totals):
                contributions = mapping.sum(axis=0) * baseline
                assigned += total * np.asarray([contributions[mask].sum() for mask in weights[:2]]) / contributions.sum()
            masses = [baseline[mask].sum() for mask in weights]
            header = dict(source=args.source, gene_id=gene_id, feature_id=event_record.feature_id, event_type=event.event_type, n_contrasts=event_record.n_contrasts, n_transcripts=len(transcripts), n_ecs=event_record.n_ecs, pooled_baseline_converged=converged, pooled_event_mass=masses[0] + masses[1], pooled_inclusion=masses[0] / (masses[0] + masses[1]), gene_molecules=primer_totals.sum(), expected_included_origins=assigned[0], expected_excluded_origins=assigned[1], minimum_expected_path_origins=min(assigned))
            anchors = [(f"interior_inclusion_{psi}", binary_event_mixture(paths, inclusion=psi), interior_totals) for psi in (.2, .5, .8)]
            anchors.append(("pooled_label_independent", baseline, pooled_totals))
            for name, mixture, totals in anchors:
                information = mixture_event_information(base.compatibility, paths, mixture, totals)
                fixed, free = information["fixed_mixture_information"], information["free_transcript_information"]
                ratio = free / fixed if fixed > np.finfo(float).tiny else np.nan
                if free < 0 or fixed < 0 or (fixed > 0 and free > fixed * (1 + 1e-7) + 1e-12):
                    raise ValueError("profiling cannot add likelihood information")
                output.append({**header, "anchor": name, **information, "retained_information_fraction": ratio})
        if gene_number % 50 == 0:
            print(f"{args.source}, shard{args.shard_index}, genes={gene_number + 1}, records={len(output)}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output_dir / "events.tsv.gz", sep="\t", index=False)
    (args.output_dir / "settings.json").write_text(json.dumps(dict(source=args.source, candidate_settings=settings, shard_count=args.shard_count, requested_features=family.feature_id.tolist(), anchors=["interior_inclusion_0.2", "interior_inclusion_0.5", "interior_inclusion_0.8", "pooled_label_independent"], interior_event_mass_if_outside=.7, interior_total_per_type=1, primer_balance="actual pooled gene primer fractions, same at both virtual types", pooled_total_per_type="half of actual pooled gene primer totals", mixture_floor=1e-12, selection="entire supported catalog from frozen coverage audit, no significance or LR selection", limitation="expected local count information, not observed subject null fits, calibrated power or LR replication; assigned origins are model-based, not uniquely informative reads", production_changes=False), indent=2) + "\n")


if __name__ == "__main__":
    main()
