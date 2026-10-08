"""Test EC categorical-mixture compatibility over the complete screened genes."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts
from tealeaf.sc.ec_glmm import subset_gene_data
from tealeaf.sc.ec_diagnostics import pooled_ec_mixture_diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary", "prepared_weighted"), required=True)
    parser.add_argument("--data-cache", type=Path, help="Required only for the existing prepared-weighted input diagnostic.")
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--alias-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=8)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid EC-model diagnostic shard")
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    if settings["subject_fold"] is not None:
        raise ValueError("use the frozen full-data gene family")
    if (args.source == "prepared_weighted") != (args.data_cache is not None):
        raise ValueError("explicit data cache is only required for prepared_weighted")
    data_cache = args.data_cache if args.data_cache is not None else args.cache / f"{args.source}_paired/prepared.pkl"
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(data_cache, settings)
    counts = tuple(value.tocsc() for value in counts)
    family = pd.read_csv(args.alias_audit, sep="\t")
    family_source = "parsimony_binary" if args.source == "prepared_weighted" else args.source
    family = family.loc[family.source.eq(family_source)]
    gene_ids = sorted(family.gene_id.unique())
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    output = []
    for gene_id in gene_ids[args.shard_index::args.shard_count]:
        gene = lookup[canonical(gene_id)]
        transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
        if len(transcripts) >= 2:
            base, _, _ = subset_gene_data(counts, designs, transcripts, gene_ecs[gene], np.ones((len(metadata), 1)), np.zeros(len(metadata)), drop_zero=False)
            blocks = zip(base.counts, base.compatibility)
        else:
            # Weighted support can leave one or no transcript. This audit is
            # not a GLMM, so retain those requested genes instead of deleting
            # them to satisfy the GLMM's two-transcript input requirement.
            blocks = []
            for observed, design in zip(counts, designs):
                mapping = design[gene_ecs[gene]][:, transcripts]
                supported = np.asarray(mapping.sum(axis=1)).ravel() > 0
                blocks.append((observed[:, gene_ecs[gene]][:, supported].toarray(), mapping[supported].toarray()))
        for primer, (values, mapping) in enumerate(blocks):
            raw_total = float(counts[primer][:, gene_ecs[gene]].sum())
            output.append(dict(source=args.source, gene_id=gene_id, primer=primer, gene_ec_molecules_before_mapping_support=raw_total, supported_molecule_fraction=float(values.sum() / raw_total) if raw_total > 0 else np.nan, **pooled_ec_mixture_diagnostics(values, mapping)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output_dir / "genes.tsv.gz", sep="\t", index=False)
    (args.output_dir / "settings.json").write_text(json.dumps(dict(source=args.source, data_cache=str(data_cache), family_source=family_source, candidate_settings=settings, requested_genes=gene_ids, shard_count=args.shard_count, max_iter=300, model="separate unrestricted transcript mixture per primer, normalized categorical EC columns", certificate="KL lower bound from convex first-order tangent, not assuming optimizer convergence", probability_bound="finite-sample method-of-types bound for integral independent categorical molecules, arbitrary observation-specific mixtures allowed", count_validation="every original observation must be nonnegative and integral before pooling", selection="whole frozen supported gene family, no significance or LR selection", interpretation="EC model compatibility only on supported count rows, not a differential-splicing test or replication endpoint; mapping-support losses recorded separately", production_changes=False), indent=2) + "\n")


if __name__ == "__main__":
    main()
