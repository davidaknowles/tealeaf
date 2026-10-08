"""Fixed-random gene panel for reference-only, exact-read EC opportunities."""

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts
from tealeaf.sc.ec_glmm import subset_gene_data
from tealeaf.sc.ec_diagnostics import pooled_ec_mixture_diagnostics
from tealeaf.sc.sequence_ec import exact_sequence_ec_kernels, sequence_ec_background_kernels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary"), required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--alias-audit", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fasta-index", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=8)
    parser.add_argument("--genes", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261007)
    parser.add_argument("--read-lengths", type=int, nargs="+", default=(31, 75, 150))
    parser.add_argument("--end-window", type=int, default=1000)
    parser.add_argument("--max-sequence-bases", type=int, default=5_000_000)
    parser.add_argument("--whole-catalog", action="store_true")
    parser.add_argument("--uniform-only", action="store_true")
    parser.add_argument("--background-fraction", type=float, default=0.)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count or args.genes <= 0 or min(args.read_lengths) <= 0 or len(set(args.read_lengths)) != len(args.read_lengths) or args.end_window <= 0 or not np.isfinite(args.background_fraction) or not 0 <= args.background_fraction < 1:
        raise ValueError("invalid opportunity-audit configuration")
    with args.candidate_cache.open("rb") as handle:
        recipe = pickle.load(handle)["settings"]
    if recipe["subject_fold"] is not None:
        raise ValueError("reference gene family uses the full-data screening recipe")
    source = args.cache / f"{args.source}_paired"
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(source / "prepared.pkl", recipe)
    counts = tuple(value.tocsc() for value in counts)
    features = (source / "features.txt").read_text().splitlines()
    family = pd.read_csv(args.alias_audit, sep="\t")
    family = family.loc[family.source.eq(args.source)].groupby("gene_id").transcripts.first()
    rng = np.random.default_rng(args.seed)
    requested = sorted(family.index) if args.whole_catalog else sorted(set(rng.choice(family.index.to_numpy(), min(args.genes, len(family)), replace=False)) | {family.idxmin(), family.idxmax()})
    gene_lookup = {canonical(value): index for index, value in enumerate(genes)}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    windows = (None,) if args.uniform_only else (None, args.end_window)
    settings = dict(source=args.source, candidate_settings=recipe, requested_genes=requested, seed=args.seed, read_lengths=args.read_lengths, terminal_start_windows=["full" if value is None else value for value in windows], shard_count=args.shard_count, max_sequence_bases=args.max_sequence_bases, background_fraction=args.background_fraction, selection="whole screened gene catalog, no significance or LR selection" if args.whole_catalog else "fixed random screened genes plus minimum/maximum transcript dimensions, no significance or LR selection", model="gene-local exact unstranded single-read matching, uniform RH starts and uniform DT terminal-window starts", background="fixed compatibility-uniform per-transcript component with unchanged DT/RH sequence capture exposures, not an estimated sequencing-error rate", limitations="not the real aligner or compound-UMI EC construction; omits other-gene mapping, empirical primer bias and genomic variants; shorter seed-length control is not full-read alignment", inference="opportunity/support and fixed-map KL diagnostics only, not a DS test, full-family benchmark or replication improvement", production_changes=False)
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    rows = []
    with Fasta(str(args.fasta), indexname=str(args.fasta_index), as_raw=True, build_index=False, rebuild=False) as reference:
        names = {}
        for name in reference.keys():
            key = canonical(name)
            if key in names:
                raise ValueError("ambiguous canonical FASTA transcript identifier")
            names[key] = name
        for gene_id in requested[args.shard_index::args.shard_count]:
            gene = gene_lookup[canonical(gene_id)]
            transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
            base, _, _ = subset_gene_data(counts, designs, transcripts, gene_ecs[gene], np.ones((len(metadata), 1)), np.zeros(len(metadata)), drop_zero=False)
            if not np.array_equal(base.compatibility[0] > 0, base.compatibility[1] > 0):
                raise ValueError("binary primer supports are not identical")
            membership = base.compatibility[0] > 0
            sequences = None
            sequence_error = None
            try:
                sequences = [reference[names[canonical(features[index])]][:] for index in transcripts]
                if sum(map(len, sequences)) > args.max_sequence_bases:
                    raise ValueError("declared per-gene sequence resource cap exceeded")
            except (KeyError, ValueError) as exception:
                sequence_error = str(exception)
            for read_length in args.read_lengths:
                for window in windows:
                    started = time.perf_counter()
                    error = sequence_error
                    maps = matched = None
                    try:
                        if error is None:
                            kernel = exact_sequence_ec_kernels(sequences, read_length, end_window=window)
                            maps, matched = sequence_ec_background_kernels(kernel, membership, args.background_fraction)
                            suffix = "" if args.background_fraction == 0 else f"_E{args.background_fraction}"
                            np.savez_compressed(args.output_dir / f"{gene_id}_L{read_length}_W{window or 'full'}{suffix}.npz", dt=maps[0], rh=maps[1], matched_classes=matched, transcripts=np.asarray([features[index] for index in transcripts]), ec_indices=gene_ecs[gene], dt_valid_starts=kernel.dt_starts, rh_valid_starts=kernel.rh_starts, background_fraction=args.background_fraction)
                    except ValueError as exception:
                        error = str(exception)
                    for primer, observed in enumerate(base.counts):
                        total = float(observed.sum())
                        row = dict(source=args.source, gene_id=gene_id, primer=primer, read_length=read_length, end_window="full" if window is None else str(window), background_fraction=args.background_fraction, transcripts=len(transcripts), ecs=observed.shape[1], input_molecules=total, elapsed_kernel_seconds=time.perf_counter() - started, error=error or "", status="sequence_or_class_failure" if error else "ok")
                        if error is None:
                            mapping = maps[primer]
                            supported = mapping.sum(axis=1) > 0
                            represented = float(observed[:, supported].sum())
                            row.update(exact_class_fraction=float(matched.mean()), retained_molecule_fraction=represented / total if total > 0 else np.nan, unrepresented_molecules=total - represented, full_positive_support=represented == total and total > 0, generated_valid_starts=int((kernel.dt_starts if primer == 0 else kernel.rh_starts).sum()))
                            if represented > 0:
                                # Both discrepancies use exactly the same supported
                                # rows. Missing molecules remain explicit and must
                                # not be credited as a full-model goodness-of-fit.
                                exact = pooled_ec_mixture_diagnostics(observed[:, supported], mapping[supported])
                                old = pooled_ec_mixture_diagnostics(observed[:, supported], base.compatibility[primer][supported])
                                row.update(sequence_KL_lower_bound=exact["KL_lower_bound"], baseline_same_rows_KL_lower_bound=old["KL_lower_bound"], sequence_converged=exact["fit_converged"], baseline_same_rows_converged=old["fit_converged"], sequence_KL=exact["KL"], baseline_same_rows_KL=old["KL"])
                        rows.append(row)
                    print(f"{args.source}, {gene_id}, L={read_length}, W={window}, status={error or 'ok'}", flush=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "genes.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
