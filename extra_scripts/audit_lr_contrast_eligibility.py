#!/usr/bin/env python3
"""Audit exact pre-fit eligibility of frozen native LR rank-prefix contrasts."""

import argparse
import json
from pathlib import Path

import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts, event_path_index
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, modeled_gene_umis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published-ranks", type=Path, required=True)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--coverage-manifest", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new audit output required")
    ranked = pd.read_csv(args.published_ranks, sep="\t", low_memory=False)
    ranked = ranked.loc[ranked.method.eq("SUPPA2 (full data)")].copy()
    if len(ranked) != 200 or set(ranked['rank']) != set(range(1, 201)):
        raise ValueError("complete unchanged native top-200 prefix required")
    settings = json.loads(args.coverage_manifest.read_text())["candidate_settings"]
    if settings["subject_fold"] is not None:
        raise ValueError("original full-data coverage recipe required")
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(args.control / "prepared.pkl", settings)
    counts = tuple(value.tocsc() for value in counts)
    features = (args.control / "features.txt").read_text().splitlines()
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    contexts, records = {}, []
    for row in ranked.itertuples(index=False):
        event = catalog.loc[row.feature_id]
        gene_id = canonical(event.gene_id)
        if gene_id not in contexts:
            gene = lookup.get(gene_id)
            status, transcripts, pairs = "gene absent from prepared inputs", (), {}
            if gene is not None:
                if not 0 < len(gene_ecs[gene]) <= settings["max_ecs"]:
                    status = "gene EC-count screen"
                else:
                    transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
                    status = "gene transcript-support screen"
                    if len(transcripts) >= 2:
                        umis = modeled_gene_umis(counts, designs, gene_ecs[gene], transcripts)
                        specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=settings["min_gene_umis"], min_samples=settings["min_gene_samples"], min_celltype_mice=settings["min_celltype_mice"])
                        pairs = {tuple(sorted(map(str, coverage[-1]))): int(metadata.iloc[coverage[0]].mouse.nunique()) for coverage, _ in specs}
                        status = "contrast gene-coverage or subject-count screen"
            contexts[gene_id] = status, transcripts, pairs
        status, transcripts, pairs = contexts[gene_id]
        key = tuple(sorted((str(row.level_a), str(row.level_b))))
        subjects = pairs.get(key, 0)
        if key in pairs:
            status = "requested before fitting" if event_path_index(transcripts, features, event.included, event.excluded) is not None else "event class-support screen"
        records.append(dict(rank=row.rank, feature_id=row.feature_id, contrast_id=row.contrast_id, gene_id=gene_id, eligibility=status, requested_subjects=subjects))
    table = pd.DataFrame(records)
    summaries = [dict(cutoff=cutoff, eligibility=label, requested=len(local), unique_genes=local.gene_id.nunique(), median_requested_subjects=float(local.requested_subjects.median())) for cutoff in (100, 200) for label, local in table.loc[table['rank'].le(cutoff)].groupby("eligibility")]
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "rank_prefix.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    inputs = (args.published_ranks, args.coverage_manifest, args.event_catalog, args.control / "manifest.json", args.control / "features.txt")
    manifest = dict(input_hashes={str(path): file_hash(path) for path in inputs}, candidate_settings=settings, scope="exact pre-fit screening of unchanged published native LR prefixes; no new ranking, fitting, event selection, discovery or replication claim", limitation="Requested does not imply numerical fit success, valid p-values or an available independent report. The old full sequence fit assessment remains pending.", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
