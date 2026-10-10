#!/usr/bin/env python3
"""Binary-path EC random-subject likelihood-ratio tests for pairwise contrasts.

Reads the same prepared EC cache and candidate cache as run_paired_path_test.py
and writes the same shard layout (paired_path.tsv, path_usage.tsv,
summary.json), so the split and long-read assessment scripts apply unchanged.
Only two-path blocks are fitted; multi-path candidates are written as
unfitted (p = 1) so the declared family is preserved.
"""

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc import ec_glmm
from tealeaf.sc.differential import helmert_basis
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.junction_benchmark import JunctionBundle, benjamini_hochberg
from tealeaf.sc.junction_paths import block_path_exons, path_junction_map
from tealeaf.sc.path_marginal import BinaryECPathLikelihood
from tealeaf.sc.local_path_reads import path_read_opportunities, pooled_path_shares
from extra_scripts.collect_local_path_reads import path_key

from tealeaf.sc.local_read_mixed import LocalReadMixed, local_read_mixed_test

PRIMERS = ("poly(dT)", "random hexamer")
READ_LENGTH = 151


def local_class_masks(opportunities, n_paths, anchors):
    """Observed read classes; without anchors the all-paths class is dropped."""
    full = (1 << n_paths) - 1
    return sorted(mask for mask in opportunities if anchors or mask != full)


def local_read_likelihood(lookup, block, opportunities, subjects, labels, levels, anchors, target=0, shares=None):
    """Primer-specific block-local class counts, target path versus the rest.

    The compatibility of class k with the target is its read-opportunity count
    n[k][target]; with the collapsed rest it is sum_j w_j n[k][j] over the
    other paths, w their label-blind pooled shares renormalized. For two paths
    this is the plain binary model. Normalizers are effective local lengths.
    """
    n_paths = len(next(iter(opportunities.values())))
    masks = local_class_masks(opportunities, n_paths, anchors)
    rest = np.ones(n_paths) if shares is None else np.array(shares, dtype=float)
    rest[target] = 0.
    rest /= rest.sum()
    components = np.column_stack([np.zeros(len(masks)), [opportunities[mask][target] for mask in masks], [opportunities[mask] @ rest for mask in masks]])
    if (components[:, 1].sum() == 0) or (components[:, 2].sum() == 0):
        raise ValueError("a path has no read opportunity in the chosen classes")
    keys, rows = [], []
    for subject in np.unique(subjects):
        for label in np.unique(labels[subjects == subject]):
            values = [np.array([lookup.get((block, subject, levels[int(label)], primer, mask), 0.) for mask in masks]) for primer in PRIMERS]
            if sum(value.sum() for value in values) > 0:
                keys.append((subject, label))
                rows.append(values)
    if not rows:
        raise ValueError("no block-local molecules")
    counts = tuple(np.asarray([row[primer] for row in rows]) for primer in range(len(PRIMERS)))
    n = len(keys)
    return BinaryECPathLikelihood(counts, (components, components), np.asarray([key[0] for key in keys]), np.asarray([key[1] for key in keys]), np.ones((n, 2)), np.zeros(n))


def simes(p_values):
    """Simes combination of a block's per-path p-values."""
    ordered = np.sort(np.asarray(p_values, dtype=float))
    return float(min(1., np.min(len(ordered) * ordered / np.arange(1, len(ordered) + 1))))


def junction_likelihood(bundle, index, sample_row, block, signatures, subjects, labels, levels):
    """One 'EC' per variable junction, compatibility = path membership.

    Rows are subject-by-type junction UMI vectors; the multinomial normalizer
    sum_j (M psi)_j = |J_0| psi + |J_1| (1 - psi) is the effective-length
    correction of tealeaf.sc.junction_paths.
    """
    anchor = lambda value: None if value in (None, "null") or (isinstance(value, str) and value.strip() == "null") else tuple(json.loads(value) if isinstance(value, str) else value)
    paths = [block_path_exons(anchor(block["left_anchor"]), signature, anchor(block["right_anchor"])) for signature in signatures]
    keys, membership = path_junction_map(paths, block["strand"])
    if len(keys) == 0 or (membership.sum(axis=0) == 0).any():
        raise ValueError("a path has no variable junction")
    columns = [index.get((block["chromosome"], start + 1, end), []) for start, end in keys]
    keys_out, counts = [], []
    for subject in np.unique(subjects):
        for label in np.unique(labels[subjects == subject]):
            position = sample_row.get((subject, levels[int(label)]))
            if position is None:
                continue
            values = np.array([bundle.counts[position, column].sum() if column else 0. for column in columns])
            if values.sum() > 0:
                keys_out.append((subject, label))
                counts.append(values)
    if not counts:
        raise ValueError("no junction UMIs")
    components = np.column_stack([np.zeros(len(keys)), membership[:, 0], membership[:, 1]])
    n = len(keys_out)
    return BinaryECPathLikelihood((np.asarray(counts),), (components,), np.asarray([key[0] for key in keys_out]), np.asarray([key[1] for key in keys_out]), np.ones((n, 2)), np.zeros(n))
from tealeaf.sc.path_marginal_grid import binary_grid_test, estimate_primer_offset, prepare_grid_likelihood


def merge(root, shard_count):
    """Concatenate one cohort's shards and add BH FDR over the whole family."""
    paths = [root / f"shard_{index}" / "paired_path.tsv" for index in range(shard_count)]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"incomplete cohort, missing {missing[:3]}")
    table = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    table["fdr"] = benjamini_hochberg(table.p_value.to_numpy())
    (root / "merged").mkdir(exist_ok=True)
    table.to_csv(root / "merged" / "paired_path.tsv", sep="\t", index=False, na_rep="NA")
    print(f"{root}: {len(table)} tests, {int(table.converged.sum())} converged, {int((table.fdr < .05).sum())} BH", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--test-id-file", type=Path)
    parser.add_argument("--nodes", type=int, default=9)
    parser.add_argument("--spliced-ecs-only", action="store_true")
    parser.add_argument("--junction-prefix", type=Path, help="fit junction UMIs on the same block paths instead of ECs")
    parser.add_argument("--block-cache", type=Path)
    parser.add_argument("--local-reads", type=Path, help="directory from collect_local_path_reads.py --collate")
    parser.add_argument("--no-anchors", action="store_true", help="with --local-reads, drop the both-paths class")
    parser.add_argument("--binomial-mixed", action="store_true", help="with --local-reads, test path-0-only vs path-1-only counts with the local-read binomial random intercept/slope model")
    parser.add_argument("--multi-path", action="store_true", help="with --local-reads, fit S one-versus-rest tests per multi-path block and combine by Simes")
    parser.add_argument("--primer-offset", action="store_true", help="estimate a label-blind second-primer logit offset per test")
    parser.add_argument("--merge", action="store_true", help="merge completed shards under --output-dir instead of fitting")
    args = parser.parse_args()
    if args.merge:
        merge(args.output_dir, args.shard_count)
        return
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached["settings"]
    if settings["test_effect"] != "cell_type_pairwise":
        raise ValueError("pairwise cell-type candidates required")
    candidates = cached["candidates"]
    if args.test_id_file is not None:
        requested = set(pd.read_csv(args.test_id_file, sep="\t").test_id)
        candidates = [candidate for candidate in candidates if candidate[0] in requested]
    if args.max_candidates is not None:
        candidates = candidates[:args.max_candidates]
    candidates = partition_candidates(candidates, args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, settings)
    if args.spliced_ecs_only:
        precursor = np.array([name.endswith("-I") for name in Path(settings["features"]).read_text().splitlines()])
    if args.junction_prefix is not None:
        import gzip
        bundle = JunctionBundle.load(args.junction_prefix)
        bundle.counts = bundle.counts.tocsr()
        junction_index = {}
        for column, row in enumerate(bundle.junctions[["chromosome", "start", "end"]].itertuples(index=False)):
            junction_index.setdefault((row.chromosome, int(row.start), int(row.end)), []).append(column)
        sample_row = {(row.subject, row.cell_type): position for position, row in enumerate(bundle.samples.itertuples(index=False))}
        blocks = {row["block_id"]: row for row in json.load(gzip.open(args.block_cache or settings["block_cache"], "rt"))}
    if args.local_reads is not None:
        import gzip
        local_counts = pd.read_csv(args.local_reads / "counts.tsv.gz", sep="\t")
        local_lookup = {(row.path_key, row.subject, row.cell_type, row.primer, int(row.mask)): float(row.count) for row in local_counts.itertuples(index=False)}
        local_blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
        local_opportunities = {}
    rows_out, usage_out, baselines = [], [], {}
    started = time.monotonic()
    for candidate in candidates:
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, levels = candidate
        record = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "contrast": "cell_type_pairwise", "level_a": levels[0], "level_b": levels[1], "method": "binary_ec_random_subject_grid", "n_paths": len(signatures), "path_signatures": json.dumps(signatures), "n_subjects": 0, "median_gene_umis": np.nan, "statistic": 0., "p_value": 1., "raw_p_value": 1., "converged": False, "mean_difference": "[]", "mean_difference_norm": np.nan, "error": ""}
        begin = time.monotonic()
        try:
            if len(signatures) != 2 and not args.multi_path:
                raise ValueError("multi-path block, not fitted by the binary model")
            local_metadata, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
            subjects = local_metadata.mouse.astype(str).to_numpy()
            base, _, totals = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
            record["median_gene_umis"] = float(np.median(totals))
            local_index = np.asarray(path_index)
            if args.binomial_mixed:
                block_key = path_key(block_id, signatures)
                tensor = np.asarray([[[[local_lookup.get((block_key, subject, levels[c], primer, mask), 0.) for mask in (1, 2)] for c in (0, 1)] for primer in PRIMERS] for subject in np.unique(subjects)])
                mixed = local_read_mixed_test(LocalReadMixed(tensor))
                effect = mixed["log_odds_effect"]
                record.update({"n_subjects": mixed["n_subjects"], "statistic": mixed["statistic"], "p_value": mixed["p_value"], "raw_p_value": mixed["p_value"], "converged": mixed["converged"], "mean_difference": json.dumps([effect / np.sqrt(2)]), "mean_difference_norm": abs(effect) / np.sqrt(2)})
                if mixed["converged"]:
                    for c, sign in ((0, -.5), (1, .5)):
                        proportion = float(1 / (1 + np.exp(-sign * effect)))
                        for path, value in enumerate((proportion, 1 - proportion)):
                            usage_out.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "subject": "model", "cell_type": levels[c], "path": f"Path {path + 1}", "path_number": path + 1, "path_signature": json.dumps(signatures[path]), "proportion": value})
                record["runtime_seconds"] = time.monotonic() - begin
                rows_out.append(record)
                continue
            if args.multi_path and len(signatures) > 2:
                block_key = path_key(block_id, signatures)
                if block_key not in local_opportunities:
                    local_opportunities[block_key] = path_read_opportunities(local_blocks[block_key]["paths"], READ_LENGTH)
                opportunities = local_opportunities[block_key]
                masks = local_class_masks(opportunities, len(signatures), not args.no_anchors)
                pooled = {mask: sum(local_lookup.get((block_key, subject, levels[int(label)], primer, mask), 0.) for subject, label in set(zip(subjects, labels)) for primer in PRIMERS) for mask in masks}
                shares = pooled_path_shares(pooled, {mask: opportunities[mask] for mask in masks})
                p_values, statistics, effects, means = [], [], np.full(len(signatures), np.nan), {}
                for target in range(len(signatures)):
                    try:
                        likelihood = local_read_likelihood(local_lookup, block_key, opportunities, subjects, labels, levels, not args.no_anchors, target, shares)
                        if args.primer_offset:
                            likelihood.primer_offsets = np.array([0., estimate_primer_offset(likelihood)])
                        fit = binary_grid_test(likelihood, nodes=args.nodes)
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
                        p_values.append(1.)
                        continue
                    p_values.append(fit["p_value"])
                    statistics.append(fit["statistic"])
                    if fit["converged"]:
                        effects[target] = fit["standardized_means"][1, 0] - fit["standardized_means"][0, 0]
                        means[target] = fit["standardized_means"][:, 0]
                record.update({"n_subjects": len(np.unique(subjects)), "statistic": max(statistics, default=0.), "p_value": simes(p_values), "converged": bool(np.isfinite(effects).any()), "n_path_tests": len(signatures), "n_converged_path_tests": int(np.isfinite(effects).sum()), "path_shares": json.dumps(shares.tolist())})
                record["raw_p_value"] = record["p_value"]
                if record["converged"]:
                    centered = np.where(np.isfinite(effects), effects, 0.)
                    record["mean_difference"] = json.dumps((helmert_basis(len(signatures)).T @ centered).tolist())
                    record["mean_difference_norm"] = float(np.linalg.norm(centered))
                    for target, values in means.items():
                        for c in (0, 1):
                            usage_out.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "subject": "model", "cell_type": levels[c], "path": f"Path {target + 1}", "path_number": target + 1, "path_signature": json.dumps(signatures[target]), "proportion": float(values[c])})
                record["runtime_seconds"] = time.monotonic() - begin
                rows_out.append(record)
                continue
            if args.local_reads is not None:
                block_key = path_key(block_id, signatures)
                if block_key not in local_opportunities:
                    local_opportunities[block_key] = path_read_opportunities(local_blocks[block_key]["paths"], READ_LENGTH)
                likelihood = local_read_likelihood(local_lookup, block_key, local_opportunities[block_key], subjects, labels, levels, not args.no_anchors)
            elif args.junction_prefix is not None:
                likelihood = junction_likelihood(bundle, junction_index, sample_row, blocks[block_id], signatures, subjects, labels, levels)
            else:
                if args.spliced_ecs_only:
                    base, kept = ec_glmm.spliced_ec_data(base, precursor[np.asarray(transcripts)])
                    local_index = local_index[kept]
                key = (gene, tuple(rows), tuple(np.asarray(transcripts)), args.spliced_ecs_only)
                if key not in baselines:
                    baseline, converged = pooled_isoform_weights(base, max_iter=250, return_status=True)
                    if not converged:
                        raise ValueError("pooled transcript baseline did not converge")
                    baselines[key] = baseline
                likelihood = prepare_grid_likelihood(base, local_index, labels, subjects, baselines[key])
            if args.primer_offset:
                likelihood.primer_offsets = np.array([0., estimate_primer_offset(likelihood)])
                record["primer_offset"] = float(likelihood.primer_offsets[1])
            result = binary_grid_test(likelihood, nodes=args.nodes)
            difference = float(result["coefficients"][1]) / np.sqrt(2)
            record.update({"n_subjects": result["n_subjects"], "statistic": result["statistic"], "p_value": result["p_value"], "raw_p_value": result["p_value"], "converged": result["converged"], "mean_difference": json.dumps([difference]), "mean_difference_norm": abs(difference), "null_concentration": result["null_concentration"], "alternative_concentration": result["alternative_concentration"], "null_subject_sd": result["null_subject_sd"], "alternative_subject_sd": result["alternative_subject_sd"], "quadrature_error": result["quadrature_error"], "null_iterations": result["null_fit"].nit, "alternative_iterations": result["alternative_fit"].nit, "final_nodes": result["nodes"], "final_refine": result["refine"]})
            if result["converged"]:
                for level, means in zip(result["levels"], result["standardized_means"]):
                    for path, proportion in enumerate(means):
                        usage_out.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "subject": "model", "cell_type": levels[int(level)], "path": f"Path {path + 1}", "path_number": path + 1, "path_signature": json.dumps(signatures[path]), "proportion": proportion})
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exception:
            record["error"] = str(exception)
        record["runtime_seconds"] = time.monotonic() - begin
        rows_out.append(record)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows_out).to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(usage_out, columns=["test_id", "block_id", "gene_id", "subject", "cell_type", "path", "path_number", "path_signature", "proportion"]).to_csv(args.output_dir / "path_usage.tsv", sep="\t", index=False, na_rep="NA")
    fitted = [row for row in rows_out if row["n_paths"] == 2]
    summary = {"candidates": len(rows_out), "binary_candidates": len(fitted), "converged": int(sum(row["converged"] for row in fitted)), "elapsed_seconds": time.monotonic() - started, "candidate_settings": settings, "nodes": args.nodes, "spliced_ecs_only": args.spliced_ecs_only, "junction_counts": args.junction_prefix is not None, "local_reads": str(args.local_reads), "anchors": not args.no_anchors, "primer_offset": args.primer_offset, "multi_path": args.multi_path, "model": "binary EC Beta random-subject GLMM, grid/AGHQ integration, LRT chi-square"}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "candidate_settings"}), flush=True)


if __name__ == "__main__":
    main()
