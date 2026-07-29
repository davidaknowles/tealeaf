#!/usr/bin/env python3
"""Estimate subisoform covariance and test differential splicing in pseudobulks."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats

from tealeaf.sc import differential, glm_cv


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument("--transcript-to-gene", required=True, type=Path)
    parser.add_argument("--fit-prefix", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--barcode-groups", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--block-cache", type=Path)
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-pseudobulk-umis", type=float, default=100_000)
    parser.add_argument("--min-gene-umis", type=float, default=25)
    parser.add_argument("--min-condition-replicates", type=int, default=3)
    parser.add_argument("--max-paths", type=int, default=8)
    parser.add_argument("--max-logratio-variance", type=float, default=0.125)
    parser.add_argument("--min-path-proportion", type=float, default=0.15)
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--batch-cells", type=int, default=128)
    parser.add_argument("--fit-max-iter", type=int, default=100)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--bootstrap-blocks", type=int, default=20)
    parser.add_argument("--bootstrap-reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def emit(event, **fields):
    print(json.dumps({"event": event, **fields}), flush=True)


def load_blocks(path, gtf):
    if path is not None and path.is_file():
        with gzip.open(path, "rt") as handle:
            payload = json.load(handle)
        return [
            differential.SpliceBlock(
                gene_id=row["gene_id"],
                block_id=row["block_id"],
                chromosome=row["chromosome"],
                strand=row["strand"],
                left_anchor=(
                    tuple(row["left_anchor"])
                    if row["left_anchor"] is not None
                    else None
                ),
                right_anchor=(
                    tuple(row["right_anchor"])
                    if row["right_anchor"] is not None
                    else None
                ),
                transcripts=tuple(row["transcripts"]),
                path_index=tuple(row["path_index"]),
                path_signatures=tuple(
                    tuple(tuple(interval) for interval in signature)
                    for signature in row["path_signatures"]
                ),
            )
            for row in payload
        ]
    blocks = differential.build_splice_blocks(gtf)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as handle:
            json.dump([
                {
                    "gene_id": block.gene_id,
                    "block_id": block.block_id,
                    "chromosome": block.chromosome,
                    "strand": block.strand,
                    "left_anchor": block.left_anchor,
                    "right_anchor": block.right_anchor,
                    "transcripts": block.transcripts,
                    "path_index": block.path_index,
                    "path_signatures": block.path_signatures,
                }
                for block in blocks
            ], handle)
    return blocks


def parse_group(value):
    parts = str(value).rsplit("__", 2)
    if len(parts) != 3:
        raise ValueError(f"group does not have cluster__DX__mouse form: {value}")
    return parts


def benjamini_hochberg(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def aggregate_inputs(args, prepared):
    pairs = pd.read_csv(args.primer_pairs, sep="\t", dtype=str)
    cell_to_poly = dict(zip(pairs["cell_id"], pairs["polydt_barcode"]))
    group_table = pd.read_csv(
        args.barcode_groups,
        header=None,
        names=["barcode", "group"],
        dtype=str,
    ).drop_duplicates("barcode")
    barcode_to_group = dict(zip(group_table["barcode"], group_table["group"]))
    assigned = np.array([
        barcode_to_group.get(cell_to_poly.get(str(cell)), None)
        for cell in prepared.barcodes
    ], dtype=object)
    raw = prepared.cv_raw_counts.tocsr()
    n_ec = raw.shape[1] // 2
    cell_totals = np.asarray(raw.sum(axis=1)).ravel()
    candidate_groups = sorted({value for value in assigned if value is not None})
    candidate_index = {value: index for index, value in enumerate(candidate_groups)}
    initial = np.array(
        [candidate_index.get(value, -1) for value in assigned], dtype=np.int64
    )
    cell_counts = np.bincount(
        initial[initial >= 0], minlength=len(candidate_groups)
    )
    umi_counts = np.bincount(
        initial[initial >= 0],
        weights=cell_totals[initial >= 0],
        minlength=len(candidate_groups),
    )
    keep = (
        (cell_counts >= int(args.min_cells))
        & (umi_counts >= float(args.min_pseudobulk_umis))
    )
    groups = [group for group, retained in zip(candidate_groups, keep) if retained]
    old_to_new = np.full(len(candidate_groups), -1, dtype=np.int64)
    old_to_new[np.flatnonzero(keep)] = np.arange(keep.sum())
    group_index = np.full(len(initial), -1, dtype=np.int64)
    valid = initial >= 0
    group_index[valid] = old_to_new[initial[valid]]
    valid = group_index >= 0
    selector = sp.coo_matrix(
        (
            np.ones(valid.sum(), dtype=np.float64),
            (group_index[valid], np.flatnonzero(valid)),
        ),
        shape=(len(groups), len(group_index)),
    ).tocsr()
    counts = (
        (selector @ raw[:, :n_ec]).tocsr(),
        (selector @ raw[:, n_ec:]).tocsr(),
    )
    return groups, group_index, cell_totals, counts


def load_aggregated_theta(args, prepared, groups, group_index, cell_totals):
    fit_rows = np.loadtxt(f"{args.fit_prefix}glm_rows.txt", dtype=str)
    fit_transcripts = np.loadtxt(f"{args.fit_prefix}glm_cols.txt", dtype=str)
    if not np.array_equal(fit_transcripts, np.asarray(prepared.features, dtype=str)):
        raise ValueError("hierarchical fit transcripts do not match prepared data")
    prepared_position = {
        str(cell): index for index, cell in enumerate(prepared.barcodes)
    }
    positions = np.array(
        [prepared_position.get(str(cell), -1) for cell in fit_rows], dtype=np.int64
    )
    if np.any(positions < 0):
        raise ValueError("hierarchical fit contains cells absent from preparation")
    with np.load(f"{args.fit_prefix}glm_factors.npz") as saved:
        state = saved["right"]
    with np.load(f"{args.fit_prefix}hierarchical_parameters.npz") as saved:
        parameters = {name: saved[name] for name in saved.files}
    aggregate, _ = differential.aggregate_hierarchical_theta(
        state,
        parameters,
        group_index[positions],
        cell_totals[positions],
        batch_cells=args.batch_cells,
    )
    if len(aggregate) != len(groups):
        raise ValueError("decoded group count does not match pseudobulk groups")
    return aggregate, fit_transcripts, parameters


def gene_structures(prepared, transcript_to_gene):
    transcript_gene_matrix, genes = glm_cv._transcript_gene_assignment(
        prepared.features, transcript_to_gene
    )
    transcript_gene = np.asarray(
        transcript_gene_matrix.indices, dtype=np.int64
    )
    n_ec = prepared.cv_raw_counts.shape[1] // 2
    designs = tuple(
        (2.0 * prepared.compatibility[start : start + n_ec]).tocsr()
        for start in (0, n_ec)
    )
    support = designs[0].copy()
    support.data.fill(1.0)
    other = designs[1].copy()
    other.data.fill(1.0)
    support = (support + other).tocsr()
    support.data.fill(1.0)
    ec_gene = glm_cv._unambiguous_ec_gene_assignment(
        support, transcript_gene_matrix
    ).tocsc()
    gene_transcripts = [
        np.flatnonzero(transcript_gene == gene) for gene in range(len(genes))
    ]
    gene_ecs = [
        ec_gene.indices[ec_gene.indptr[gene] : ec_gene.indptr[gene + 1]]
        for gene in range(len(genes))
    ]
    return np.asarray(genes, dtype=str), gene_transcripts, gene_ecs, designs


def block_mapping(block, transcript_names):
    transcript_position = {
        transcript: index for index, transcript in enumerate(transcript_names)
    }
    represented = [
        (transcript_position[transcript], path)
        for transcript, path in zip(block.transcripts, block.path_index)
        if transcript in transcript_position
    ]
    paths = sorted({path for _, path in represented})
    if len(paths) < 2:
        return None
    path_remap = {path: index for index, path in enumerate(paths)}
    mapping = np.full(len(transcript_names), -1, dtype=np.int64)
    for transcript, path in represented:
        mapping[transcript] = path_remap[path]
    signatures = [block.path_signatures[path] for path in paths]
    return mapping, signatures


def estimate_blocks(
    args,
    blocks,
    groups,
    theta,
    counts,
    genes,
    gene_transcripts,
    gene_ecs,
    designs,
):
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    feature_names = np.asarray(
        np.loadtxt(f"{args.fit_prefix}glm_cols.txt", dtype=str)
    )
    group_metadata = [parse_group(group) for group in groups]
    differential_inputs = defaultdict(list)
    bootstrap_candidates = []
    estimates_path = args.output_dir / "subisoform_estimates.jsonl.gz"
    blocks_path = args.output_dir / "subisoform_blocks.tsv"
    block_rows = []
    estimated = 0
    converged = 0
    conditional_identifiable = 0
    conditional_reliable = 0
    identifiable = 0
    profile_reliable = 0
    cached_gene = None
    cached_gene_data = None
    with gzip.open(estimates_path, "wt") as output:
        for block in blocks:
            gene = gene_to_index.get(block.gene_id)
            if gene is None:
                continue
            if gene != cached_gene:
                cached_gene = gene
                cached_gene_data = None
                transcript_indices = gene_transcripts[gene]
                ec_rows = gene_ecs[gene]
                if len(transcript_indices) < 2 or not len(ec_rows):
                    continue
                local_designs = tuple(
                    design[ec_rows][:, transcript_indices].tocsr()
                    for design in designs
                )
                local_counts = tuple(
                    matrix[:, ec_rows].tocsr() for matrix in counts
                )
                gene_umis = sum(
                    np.asarray(matrix.sum(axis=1)).ravel()
                    for matrix in local_counts
                )
                eligible = np.flatnonzero(
                    gene_umis >= float(args.min_gene_umis)
                )
                if len(eligible) < 2:
                    continue
                cached_gene_data = (
                    transcript_indices,
                    local_designs,
                    local_counts,
                    gene_umis,
                    eligible,
                )
            if cached_gene_data is None:
                continue
            (
                transcript_indices,
                local_designs,
                local_counts,
                gene_umis,
                eligible,
            ) = cached_gene_data
            mapping_result = block_mapping(block, feature_names[transcript_indices])
            if mapping_result is None:
                continue
            path_index, signatures = mapping_result
            n_paths = int(path_index.max()) + 1
            if n_paths > int(args.max_paths):
                continue
            block_rows.append({
                "block_id": block.block_id,
                "gene_id": block.gene_id,
                "chromosome": block.chromosome,
                "strand": block.strand,
                "n_paths": n_paths,
                "path_signatures": json.dumps(signatures),
                "n_eligible_pseudobulks": len(eligible),
                "total_gene_umis": float(gene_umis.sum()),
            })
            for sample in eligible:
                observed = tuple(
                    matrix.getrow(sample).toarray().ravel()
                    for matrix in local_counts
                )
                baseline = np.maximum(
                    theta[sample, transcript_indices].astype(float), 1e-12
                )
                fit = differential.fit_path_perturbation(
                    observed,
                    local_designs,
                    baseline,
                    path_index,
                    max_iter=args.fit_max_iter,
                )
                profile = differential.profiled_path_covariance(
                    fit.theta,
                    path_index,
                    local_designs,
                    tuple(value.sum() for value in observed),
                )
                estimated += 1
                converged += int(fit.converged)
                conditional_identifiable += int(fit.covariance.identifiable)
                identifiable += int(profile.identifiable)
                conditional_ok = bool(
                    fit.covariance.identifiable
                    and np.min(fit.path_proportions)
                    >= args.min_path_proportion
                    and np.linalg.eigvalsh(
                        fit.covariance.covariance
                    ).max() <= args.max_logratio_variance
                )
                profile_ok = bool(
                    profile.identifiable
                    and np.min(fit.path_proportions)
                    >= args.min_path_proportion
                    and np.linalg.eigvalsh(profile.covariance).max()
                    <= args.max_logratio_variance
                )
                conditional_reliable += int(conditional_ok)
                profile_reliable += int(profile_ok)
                cluster, condition, mouse = group_metadata[sample]
                record = {
                    "block_id": block.block_id,
                    "gene_id": block.gene_id,
                    "group": groups[sample],
                    "cluster": cluster,
                    "condition": condition,
                    "mouse": mouse,
                    "gene_umis": float(gene_umis[sample]),
                    "path_proportions": fit.path_proportions.tolist(),
                    "path_logratios": fit.path_logratios.tolist(),
                    "conditional_covariance": fit.covariance.covariance.tolist(),
                    "conditional_identifiable": fit.covariance.identifiable,
                    "conditional_information_rank": (
                        fit.covariance.information_rank
                    ),
                    "conditional_reliable": conditional_ok,
                    "profile_covariance": profile.covariance.tolist(),
                    "profile_identifiable": profile.identifiable,
                    "profile_reliable": profile_ok,
                    "profile_information_rank": profile.information_rank,
                    "profile_null_projection": profile.null_projection,
                    "fit_converged": fit.converged,
                    "fit_iterations": fit.iterations,
                }
                output.write(json.dumps(record) + "\n")
                if fit.converged and conditional_ok:
                    differential_inputs[
                        (block.block_id, cluster, "conditional")
                    ].append(
                        (
                            fit.path_logratios,
                            fit.covariance.covariance,
                            condition,
                        )
                    )
                    if profile_ok:
                        differential_inputs[
                            (block.block_id, cluster, "profile")
                        ].append(
                            (
                                fit.path_logratios,
                                profile.covariance,
                                condition,
                            )
                        )
                if (
                    fit.converged
                    and conditional_ok
                    and len(bootstrap_candidates) < args.bootstrap_blocks
                ):
                    bootstrap_candidates.append(
                        (
                            block.block_id,
                            groups[sample],
                            fit,
                            observed,
                            local_designs,
                            path_index,
                        )
                    )
            emit(
                "block_complete",
                block_id=block.block_id,
                eligible_pseudobulks=len(eligible),
                estimates=estimated,
            )
    pd.DataFrame(block_rows).to_csv(blocks_path, sep="\t", index=False)
    return differential_inputs, bootstrap_candidates, {
        "blocks": len(block_rows),
        "estimates": estimated,
        "converged": converged,
        "converged_fraction": converged / max(estimated, 1),
        "conditional_identifiable": conditional_identifiable,
        "conditional_identifiable_fraction": (
            conditional_identifiable / max(estimated, 1)
        ),
        "conditional_reliable": conditional_reliable,
        "conditional_reliable_fraction": conditional_reliable / max(estimated, 1),
        "profile_identifiable": identifiable,
        "profile_identifiable_fraction": identifiable / max(estimated, 1),
        "profile_reliable": profile_reliable,
        "profile_reliable_fraction": profile_reliable / max(estimated, 1),
    }


def differential_tests(args, inputs, rng):
    rows = []
    null_p_values = defaultdict(list)
    for (block_id, cluster, covariance_method), records in inputs.items():
        dimensions = {len(record[0]) for record in records}
        if len(dimensions) != 1:
            continue
        conditions = sorted({record[2] for record in records})
        counts = {
            condition: sum(record[2] == condition for record in records)
            for condition in conditions
        }
        conditions = [
            condition
            for condition in conditions
            if counts[condition] >= int(args.min_condition_replicates)
        ]
        retained = [record for record in records if record[2] in conditions]
        if len(conditions) < 2:
            continue
        condition_index = {
            condition: index for index, condition in enumerate(conditions)
        }
        design = np.zeros((len(retained), len(conditions)), dtype=float)
        design[:, 0] = 1
        for row, record in enumerate(retained):
            index = condition_index[record[2]]
            if index:
                design[row, index] = 1
        values = np.asarray([record[0] for record in retained], dtype=float)
        covariances = np.asarray([record[1] for record in retained], dtype=float)
        finite = np.isfinite(values).all(axis=1)
        finite &= np.isfinite(covariances).all(axis=(1, 2))
        if not finite.all():
            values = values[finite]
            covariances = covariances[finite]
            design = design[finite]
        if len(values) <= design.shape[1]:
            continue
        result = differential.multivariate_gls_test(
            values, covariances, design, tested_columns=range(1, len(conditions))
        )
        row_result = {
            "block_id": block_id,
            "cluster": cluster,
            "covariance_method": covariance_method,
            "conditions": ",".join(conditions),
            "n_samples": len(retained),
            "logratio_dimension": values.shape[1],
            "statistic": result["statistic"],
            "degrees_of_freedom": result["degrees_of_freedom"],
            "asymptotic_p_value": result["p_value"],
            "biological_variance": result["biological_variance"],
        }
        permutation_statistics = []
        for _ in range(int(args.permutations)):
            permutation = rng.permutation(len(retained))
            permuted_design = design[permutation]
            null = differential.multivariate_gls_test(
                values,
                covariances,
                permuted_design,
                tested_columns=range(1, len(conditions)),
                biological_variance=result["biological_variance"],
            )
            permutation_statistics.append(null["statistic"])
            null_p_values[covariance_method].append(null["p_value"])
        row_result["permutation_p_value"] = (
            (1 + np.sum(np.asarray(permutation_statistics) >= result["statistic"]))
            / (len(permutation_statistics) + 1)
        )
        row_result["p_value"] = result["p_value"]
        rows.append(row_result)
    table = pd.DataFrame(rows)
    if len(table):
        table["fdr"] = np.nan
        for method, positions in table.groupby("covariance_method").groups.items():
            table.loc[positions, "fdr"] = benjamini_hochberg(
                table.loc[positions, "p_value"].to_numpy()
            )
        table = table.sort_values("p_value")
    table.to_csv(args.output_dir / "differential_splicing.tsv", sep="\t", index=False)
    summary = {
        "tests": int(len(table)),
        "fdr_0.05": int((table["fdr"] <= 0.05).sum()) if len(table) else 0,
    }
    for method, values in null_p_values.items():
        values = np.asarray(values, dtype=float)
        summary[f"{method}_permutation_tests"] = int(len(values))
        summary[f"{method}_permutation_p_lt_0.05"] = float(np.mean(values < 0.05))
        summary[f"{method}_permutation_ks_p_value"] = float(
            scipy.stats.kstest(values, "uniform").pvalue
        )
    return summary


def bootstrap_validation(args, candidates, rng):
    squared_errors = []
    predicted_variances = []
    covered = []
    converged = 0
    candidate_rows = []
    for block_id, group, fit, observed, designs, path_index in candidates:
        truth = fit.path_logratios
        variance = np.diag(fit.covariance.covariance)
        candidate_squared_errors = []
        candidate_covered = []
        totals = [int(round(value.sum())) for value in observed]
        probabilities = []
        for design in designs:
            expected = np.asarray(design @ fit.theta).ravel()
            probabilities.append(expected / expected.sum())
        for _ in range(int(args.bootstrap_reps)):
            simulated = tuple(
                rng.multinomial(total, probability)
                for total, probability in zip(totals, probabilities)
            )
            replicate = differential.fit_path_perturbation(
                simulated,
                designs,
                fit.theta,
                path_index,
                max_iter=args.fit_max_iter,
            )
            if not replicate.converged:
                continue
            converged += 1
            error = replicate.path_logratios - truth
            squared_errors.extend(np.square(error))
            predicted_variances.extend(variance)
            covered.extend(np.abs(error) <= 1.96 * np.sqrt(variance))
            candidate_squared_errors.extend(np.square(error))
            candidate_covered.extend(
                np.abs(error) <= 1.96 * np.sqrt(variance)
            )
        candidate_squared_errors = np.asarray(
            candidate_squared_errors, dtype=float
        )
        candidate_rows.append({
            "block_id": block_id,
            "group": group,
            "total_umis": sum(totals),
            "minimum_path_proportion": float(
                np.min(fit.path_proportions)
            ),
            "maximum_standard_error": float(np.sqrt(variance).max()),
            "minimum_information_eigenvalue": (
                fit.covariance.minimum_positive_eigenvalue
            ),
            "bootstrap_fits": (
                len(candidate_squared_errors) // len(variance)
            ),
            "mean_squared_z": (
                float(np.mean(candidate_squared_errors / np.tile(
                    variance,
                    len(candidate_squared_errors) // len(variance),
                )))
                if len(candidate_squared_errors)
                else None
            ),
            "coverage_95": (
                float(np.mean(candidate_covered))
                if candidate_covered
                else None
            ),
        })
    pd.DataFrame(candidate_rows).to_csv(
        args.output_dir / "bootstrap_validation.tsv",
        sep="\t",
        index=False,
    )
    squared_errors = np.asarray(squared_errors, dtype=float)
    predicted_variances = np.asarray(predicted_variances, dtype=float)
    standardized_squared_errors = (
        squared_errors / predicted_variances
        if len(squared_errors)
        else np.array([], dtype=float)
    )
    return {
        "bootstrap_blocks": len(candidates),
        "bootstrap_converged_fits": converged,
        "bootstrap_variance_ratio": (
            float(squared_errors.mean() / predicted_variances.mean())
            if len(squared_errors)
            else None
        ),
        "bootstrap_mean_squared_z": (
            float(standardized_squared_errors.mean())
            if len(standardized_squared_errors)
            else None
        ),
        "bootstrap_median_squared_z": (
            float(np.median(standardized_squared_errors))
            if len(standardized_squared_errors)
            else None
        ),
        "bootstrap_component_coverage_95": (
            float(np.mean(covered)) if covered else None
        ),
    }


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    started = time.perf_counter()
    emit("preparation_started")
    prepared = glm_cv.prepare_paired_primer_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        args.primer_pairs,
        ec_design="weighted",
        regularization_target="theta",
        min_eq=5,
        min_half_umis=args.min_half_umis,
        primer_sampling_model="oligodt_tpm",
    )
    groups, group_index, cell_totals, counts = aggregate_inputs(args, prepared)
    emit(
        "pseudobulks_complete",
        seconds=time.perf_counter() - started,
        groups=len(groups),
    )
    theta, transcripts, _ = load_aggregated_theta(
        args, prepared, groups, group_index, cell_totals
    )
    emit("theta_aggregation_complete", seconds=time.perf_counter() - started)
    genes, gene_transcripts, gene_ecs, designs = gene_structures(
        prepared, args.transcript_to_gene
    )
    blocks = load_blocks(args.block_cache, args.gtf)
    if args.max_blocks is not None:
        blocks = blocks[: int(args.max_blocks)]
    emit("splice_blocks_complete", blocks=len(blocks))
    inputs, candidates, estimate_summary = estimate_blocks(
        args,
        blocks,
        groups,
        theta,
        counts,
        genes,
        gene_transcripts,
        gene_ecs,
        designs,
    )
    test_summary = differential_tests(args, inputs, rng)
    bootstrap_summary = bootstrap_validation(args, candidates, rng)
    summary = {
        "seconds": time.perf_counter() - started,
        "groups": len(groups),
        "transcripts": len(transcripts),
        **estimate_summary,
        **test_summary,
        **bootstrap_summary,
    }
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    emit("analysis_complete", **summary)


if __name__ == "__main__":
    main()
