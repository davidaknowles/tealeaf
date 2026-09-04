#!/usr/bin/env python3
"""Compare junction-level discoveries with Tealeaf in Tilgner long reads."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy.linalg import helmert
from scipy.stats import binomtest

from extra_scripts.assess_tilgner_long_read_replication import bh_adjust, read_tilgner_matrix, stable_identifier, vector_agreement, wilson_interval
from tealeaf.sc.differential import paired_mean_test
from tealeaf.sc.junction_benchmark import JunctionBundle


TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')


def mean_composition(counts):
    """Return the equal-pseudobulk mean composition for a count matrix."""
    values = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    totals = values.sum(axis=1)
    values = values[totals > 0]
    totals = totals[totals > 0]
    if not len(values):
        return np.full(values.shape[1], np.nan)
    return (values / totals[:, None]).mean(axis=0)


def mean_composition_batch(counts, sample_indices, group_indices):
    """Return equal-pseudobulk mean compositions for equal-size groups."""
    group_indices = np.asarray(group_indices, dtype=int)
    n_groups, dimension = group_indices.shape
    values = counts[sample_indices][:, group_indices.ravel()].toarray().reshape(len(sample_indices), n_groups, dimension)
    totals = values.sum(axis=2)
    proportions = np.divide(values, totals[:, :, None], out=np.full_like(values, np.nan, dtype=float), where=totals[:, :, None] > 0)
    return np.nanmean(proportions, axis=0)


def junction_from_majiq(row):
    """Convert MAJIQ exon-boundary coordinates to an intron interval."""
    left, right = sorted((int(row.start), int(row.end)))
    if right - left <= 1:
        return None
    return str(row.seqid), left + 1, right - 1, str(row.strand)


def majiq_feature_ids(table):
    """Reconstruct the stable feature identifiers used by MAJIQ normalization."""
    columns = [name for name in ("gene_id", "seqid", "strand", "event_type", "ref_exon_start", "ref_exon_end", "start", "end", "is_intron", "other_exon_start", "other_exon_end") if name in table]
    features = table[columns].astype(str).agg(":".join, axis=1)
    duplicated = features.duplicated(keep=False)
    if duplicated.any():
        features.loc[duplicated] += ":" + features.loc[duplicated].groupby(features.loc[duplicated]).cumcount().astype(str)
    return features


def load_comparators(path, mapped_cell_types, significant_only=False):
    table = pd.read_csv(path, sep="\t", low_memory=False)
    keep = (table["effect"] == "cell_type") & table["level_a"].isin(mapped_cell_types) & table["level_b"].isin(mapped_cell_types)
    if significant_only:
        keep &= table["significant"].astype(str).str.lower().eq("true")
    return table.loc[keep].copy()


def scquint_groups(selected, bundle):
    groups = {}
    metadata = bundle.junctions
    by_group = metadata[metadata["intron_group"].notna()].groupby(metadata.loc[metadata["intron_group"].notna(), "intron_group"].astype(str)).groups
    for row in selected[selected["method"] == "scQuint"].itertuples(index=False):
        indices = np.asarray(by_group.get(str(row.feature_id), []), dtype=int)
        if len(indices) < 2:
            continue
        genes = metadata.iloc[indices]["gene_id"].dropna().map(stable_identifier).unique()
        groups[("scQuint", row.contrast_id, row.feature_id)] = {"method": "scQuint", "contrast_id": row.contrast_id, "level_a": row.level_a, "level_b": row.level_b, "feature_id": row.feature_id, "p_value": row.p_value, "q_value": row.q_value, "gene_id": genes[0] if len(genes) == 1 else stable_identifier(row.gene_id), "indices": indices}
    return groups


def leafcutter_groups(selected, bundle, counts_path):
    wanted = set(selected.loc[selected["method"] == "LeafCutter", "feature_id"].astype(str))
    coordinates = defaultdict(list)
    with gzip.open(counts_path, "rt") as handle:
        next(handle)
        for line in handle:
            identifier = line.split(" ", 1)[0]
            chromosome, start, end, cluster = identifier.split(":")
            feature = f"{chromosome}:{cluster}"
            if feature in wanted:
                coordinates[feature].append((chromosome, int(start), int(end), cluster[-1]))
    lookup = {(str(row.chromosome), int(row.start), int(row.end), str(row.strand)): index for index, row in enumerate(bundle.junctions.itertuples(index=False))}
    groups = {}
    for row in selected[selected["method"] == "LeafCutter"].itertuples(index=False):
        indices = np.asarray([lookup[coordinate] for coordinate in coordinates.get(str(row.feature_id), []) if coordinate in lookup], dtype=int)
        if len(indices) < 2:
            continue
        genes = bundle.junctions.iloc[indices]["gene_id"].dropna().map(stable_identifier).unique()
        groups[("LeafCutter", row.contrast_id, row.feature_id)] = {"method": "LeafCutter", "contrast_id": row.contrast_id, "level_a": row.level_a, "level_b": row.level_b, "feature_id": row.feature_id, "p_value": row.p_value, "q_value": row.q_value, "gene_id": genes[0] if len(genes) == 1 else None, "indices": indices}
    return groups


def majiq_groups(selected, raw_directory, bundle):
    selected = selected[selected["method"] == "MAJIQ Heterogen"].copy()
    lookup = {(str(row.chromosome), int(row.start), int(row.end), str(row.strand)): index for index, row in enumerate(bundle.junctions.itertuples(index=False))}
    groups = {}
    for contrast_id, selected_contrast in selected.groupby("contrast_id"):
        raw = pd.read_csv(Path(raw_directory) / f"{contrast_id}.tsv", sep="\t", comment="#")
        raw["majiq_feature"] = majiq_feature_ids(raw)
        raw["lsv_key"] = raw[["gene_id", "event_type", "ref_exon_start", "ref_exon_end"]].astype(str).agg("|".join, axis=1)
        selected_features = set(selected_contrast["feature_id"].astype(str))
        selected_lsvs = set(raw.loc[raw["majiq_feature"].isin(selected_features), "lsv_key"])
        for lsv_key, local in raw[raw["lsv_key"].isin(selected_lsvs)].groupby("lsv_key", sort=False):
            coordinates = [junction_from_majiq(row) for row in local.itertuples(index=False)]
            indices = np.asarray(list(dict.fromkeys(lookup[coordinate] for coordinate in coordinates if coordinate in lookup)), dtype=int)
            if len(indices) < 2:
                continue
            level_a = selected_contrast["level_a"].iloc[0]
            level_b = selected_contrast["level_b"].iloc[0]
            a_column = f"{level_a}-raw_psi_quantile_0.500"
            b_column = f"{level_b}-raw_psi_quantile_0.500"
            effects = local[b_column].to_numpy(dtype=float) - local[a_column].to_numpy(dtype=float)
            effect_by_coordinate = {}
            for coordinate, effect in zip(coordinates, effects):
                if coordinate in lookup:
                    effect_by_coordinate[lookup[coordinate]] = effect
            original_delta = np.asarray([effect_by_coordinate.get(index, np.nan) for index in indices])
            if not np.isfinite(original_delta).all() or np.linalg.norm(original_delta) == 0:
                continue
            selected_local = selected_contrast[selected_contrast["feature_id"].astype(str).isin(local["majiq_feature"])]
            feature_id = f"{contrast_id}|{lsv_key}"
            groups[("MAJIQ Heterogen", contrast_id, feature_id)] = {"method": "MAJIQ Heterogen", "contrast_id": contrast_id, "level_a": level_a, "level_b": level_b, "feature_id": feature_id, "p_value": selected_local["p_value"].min(), "q_value": selected_local["q_value"].min(), "gene_id": stable_identifier(local["gene_id"].iloc[0]), "indices": indices, "majiq_delta": original_delta}
    return groups


def add_short_read_effects(groups, bundle):
    sample_groups = {cell_type: group.index.to_numpy(dtype=int) for cell_type, group in bundle.samples.groupby("cell_type")}
    batches = defaultdict(list)
    for key, group in groups.items():
        if "original_delta" not in group:
            batches[(group["level_a"], group["level_b"], len(group["indices"]))].append((key, group))
    for (level_a, level_b, _), batch in batches.items():
        indices = [group["indices"] for _, group in batch]
        first = mean_composition_batch(bundle.counts, sample_groups[level_a], indices)
        second = mean_composition_batch(bundle.counts, sample_groups[level_b], indices)
        for (_, group), delta in zip(batch, second - first):
            group["original_delta"] = delta
    output = {}
    for key, group in groups.items():
        if np.isfinite(group["original_delta"]).all() and np.linalg.norm(group["original_delta"]) > 0:
            output[key] = group
    return output


def paired_logratio_batch(counts, sample_indices, group_indices, minimum_subjects=8, pseudocount=0.5, jobs=1):
    """Fit equal-dimensional paired CLR tests after one sparse extraction."""
    import joblib

    group_indices = np.asarray(group_indices, dtype=int)
    n_groups = len(group_indices)
    if not n_groups:
        return []
    dimension = group_indices.shape[1]
    values = counts[sample_indices][:, group_indices.ravel()].toarray().reshape(len(sample_indices), n_groups, dimension)
    totals = values.sum(axis=2)
    valid = (totals[0::2] > 0) & (totals[1::2] > 0)
    logged = np.log(values + pseudocount)
    clr = logged - logged.mean(axis=2, keepdims=True)
    differences = clr[1::2] - clr[0::2]
    coordinates = np.einsum("sgk,kd->sgd", differences, helmert(dimension, full=False).T, optimize=True)

    def fit_group(group_index):
        fit = paired_mean_test(coordinates[valid[:, group_index], group_index])
        p_value = fit["p_value"] if fit["converged"] and fit["n_subjects"] >= minimum_subjects else 1.0
        return p_value, fit["n_subjects"]

    with joblib.parallel_backend("threading", n_jobs=jobs):
        return joblib.Parallel(n_jobs=jobs)(joblib.delayed(fit_group)(group_index) for group_index in range(n_groups))


def paired_logratio_two_category_batch(counts, sample_indices, group_indices, minimum_subjects=8, pseudocount=0.5):
    """Fit the two-category special case used by regression tests."""
    return paired_logratio_batch(counts, sample_indices, group_indices, minimum_subjects=minimum_subjects, pseudocount=pseudocount)


def add_paired_clr_discoveries(selected, bundle, minimum_subjects=8, jobs=1):
    candidates = selected[selected["method"] == "scQuint"].copy()
    metadata = bundle.junctions
    by_group = metadata[metadata["intron_group"].notna()].groupby(metadata.loc[metadata["intron_group"].notna(), "intron_group"].astype(str)).groups
    sample_lookup = bundle.samples.reset_index().set_index(["cell_type", "subject"])["index"]
    rows = []
    for contrast_id, local in candidates.groupby("contrast_id"):
        level_a, level_b = local[["level_a", "level_b"]].iloc[0]
        subjects = sorted(set(sample_lookup.loc[level_a].index) & set(sample_lookup.loc[level_b].index))
        sample_indices = np.asarray([sample_lookup[(level, subject)] for subject in subjects for level in (level_a, level_b)], dtype=int)
        feature_ids = local["feature_id"].astype(str).unique()
        groups_by_dimension = defaultdict(list)
        for feature_id in feature_ids:
            indices = np.asarray(by_group.get(feature_id, []), dtype=int)
            if len(indices) >= 2:
                groups_by_dimension[len(indices)].append((feature_id, indices))
        print(f"Paired CLR {contrast_id}, fitting {sum(map(len, groups_by_dimension.values())):,} groups in {len(groups_by_dimension)} dimension batches", file=sys.stderr, flush=True)
        fitted = []
        for dimension_groups in groups_by_dimension.values():
            batch = paired_logratio_batch(bundle.counts, sample_indices, [indices for _, indices in dimension_groups], minimum_subjects=minimum_subjects, jobs=jobs)
            fitted.extend((feature_id, p_value, n_subjects) for (feature_id, _), (p_value, n_subjects) in zip(dimension_groups, batch))
        if not fitted:
            continue
        frame = pd.DataFrame(fitted, columns=["feature_id", "p_value", "n_subjects"])
        frame["q_value"] = bh_adjust(frame["p_value"])
        for row in frame[frame["q_value"] < 0.05].itertuples(index=False):
            genes = metadata.iloc[np.asarray(by_group[row.feature_id], dtype=int)]["gene_id"].dropna().map(stable_identifier).unique()
            rows.append({"method": "Paired junction CLR", "contrast_id": contrast_id, "effect": "cell_type", "level_a": level_a, "level_b": level_b, "feature_id": row.feature_id, "p_value": row.p_value, "q_value": row.q_value, "significant": True, "gene_id": genes[0] if len(genes) == 1 else None})
    paired = pd.DataFrame(rows)
    return paired, scquint_groups(paired.assign(method="scQuint"), bundle)


def transcript_junction_rows(gtf_path, features, wanted):
    feature_rows = defaultdict(list)
    for row in features[features["transcript_id"].notna()].itertuples(index=False):
        feature_rows[stable_identifier(row.transcript_id)].append(int(row.row))
    result = defaultdict(list)
    last_exon = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = TRANSCRIPT_ID.search(fields[8])
            if match is None:
                continue
            transcript = stable_identifier(match.group(1))
            key = (fields[0], fields[6], transcript)
            start, end = int(fields[3]), int(fields[4])
            if key in last_exon:
                junction = (fields[0], last_exon[key] + 1, start - 1, fields[6])
                if junction in wanted:
                    result[junction].extend(feature_rows.get(transcript, []))
            last_exon[key] = end
    return {junction: np.asarray(sorted(set(rows)), dtype=int) for junction, rows in result.items()}


def source_junction_counts(groups, matrix, features, columns, gtf_path):
    wanted = {coordinate for group in groups.values() for coordinate in group["coordinates"]}
    print(f"Mapping {len(wanted):,} selected junctions to annotated source transcripts", file=sys.stderr, flush=True)
    rows = transcript_junction_rows(gtf_path, features, wanted)
    column_groups = {(cell_type, int(replicate)): group["column"].to_numpy(dtype=int) for (cell_type, replicate), group in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby(["tealeaf_cell_type", "replicate"])}
    feature_counts = {key: np.asarray(matrix[:, selected_columns].sum(axis=1)).ravel() for key, selected_columns in column_groups.items()}
    counts = {}
    for junction, feature_indices in rows.items():
        counts[junction] = {key: float(values[feature_indices].sum()) for key, values in feature_counts.items()}
    print(f"Mapped {len(rows):,} selected junctions with source-transcript support", file=sys.stderr, flush=True)
    return counts


def assess_groups(groups, source_counts):
    output = []
    for group in groups.values():
        local_coordinates = group["coordinates"]
        source = {}
        for level in (group["level_a"], group["level_b"]):
            for replicate in (1, 2):
                source[(level, replicate)] = np.asarray([source_counts.get(junction, {}).get((level, replicate), 0.0) for junction in local_coordinates])
        pooled_a = source[(group["level_a"], 1)] + source[(group["level_a"], 2)]
        pooled_b = source[(group["level_b"], 1)] + source[(group["level_b"], 2)]
        pooled_delta = pooled_b / pooled_b.sum() - pooled_a / pooled_a.sum() if pooled_a.sum() and pooled_b.sum() else np.full(len(pooled_a), np.nan)
        pooled_dot, pooled_cosine = vector_agreement(group["original_delta"], pooled_delta)
        replicate_dots = []
        for replicate in (1, 2):
            first, second = source[(group["level_a"], replicate)], source[(group["level_b"], replicate)]
            delta = second / second.sum() - first / first.sum() if first.sum() and second.sum() else np.full(len(first), np.nan)
            replicate_dots.append(vector_agreement(group["original_delta"], delta)[0])
        output.append({"method": group["method"], "contrast_id": group["contrast_id"], "level_a": group["level_a"], "level_b": group["level_b"], "feature_id": group["feature_id"], "gene_id": group["gene_id"], "p_value": group["p_value"], "q_value": group["q_value"], "n_junctions": len(local_coordinates), "mapping_complete": all(junction in source_counts for junction in local_coordinates), "short_read_effect_norm": float(np.linalg.norm(group["original_delta"])), "native_effect_dot_product": float(np.dot(group["majiq_delta"], group["original_delta"])) if "majiq_delta" in group else np.nan, "long_read_effect_norm": float(np.linalg.norm(pooled_delta)) if np.isfinite(pooled_delta).all() else np.nan, "pooled_dot_product": pooled_dot, "pooled_cosine": pooled_cosine, "pooled_replicated": pooled_dot > 0 if np.isfinite(pooled_dot) else np.nan, "replicate_1_dot_product": replicate_dots[0], "replicate_2_dot_product": replicate_dots[1], "both_replicates_replicated": replicate_dots[0] > 0 and replicate_dots[1] > 0 if np.isfinite(replicate_dots).all() else np.nan, "minimum_pooled_depth": min(pooled_a.sum(), pooled_b.sum()), "minimum_replicate_depth": min(source[(level, replicate)].sum() for level in (group["level_a"], group["level_b"]) for replicate in (1, 2))})
    return pd.DataFrame(output)


def summarize(results, minimum_depth=20, one_per_gene_contrast=False):
    rows = []
    for method, method_results in results.groupby("method"):
        complete = method_results["mapping_complete"].fillna(False).astype(bool) if "mapping_complete" in method_results else pd.Series(True, index=method_results.index)
        pooled = method_results[complete & (method_results["minimum_pooled_depth"] >= minimum_depth) & method_results["pooled_replicated"].notna()]
        strict = pooled[(pooled["minimum_replicate_depth"] >= minimum_depth / 2) & pooled["both_replicates_replicated"].notna()]
        for endpoint, group, column in (("pooled direction", pooled, "pooled_replicated"), ("both biological replicates", strict, "both_replicates_replicated")):
            if one_per_gene_contrast:
                group = group[group["gene_id"].notna()].sort_values("short_read_effect_norm", ascending=False).drop_duplicates(["gene_id", "level_a", "level_b"])
            successes = int(group[column].sum())
            if endpoint == "pooled direction":
                trials, null_rate = len(group), 0.5
            else:
                concordant = (group["replicate_1_dot_product"] > 0) == (group["replicate_2_dot_product"] > 0)
                trials = int(concordant.sum())
                null_rate = trials / (2 * len(group)) if len(group) else np.nan
            low, high = wilson_interval(successes, len(group))
            rows.append({"method": method, "endpoint": endpoint, "scope": "one call per gene and contrast" if one_per_gene_contrast else "all selected calls", "minimum_depth": minimum_depth, "n_tests": len(group), "n_replicated": successes, "replication_rate": successes / len(group) if len(group) else np.nan, "ci_low": low, "ci_high": high, "conditional_null_rate": null_rate, "orientation_p_value": binomtest(successes, trials, 0.5, alternative="greater").pvalue if trials else np.nan})
    return pd.DataFrame(rows)


def load_tealeaf_results(path):
    """Standardize the Tealeaf path-level validation table for comparison."""
    table = pd.read_csv(path, sep="\t")
    table["method"] = "Tealeaf"
    table["feature_id"] = table["test_id"]
    table["short_read_effect_norm"] = table["original_effect_norm"]
    table = table[table["mapping_complete"].fillna(False).astype(bool)].copy()
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparator-tests", required=True, type=Path)
    parser.add_argument("--junction-bundle", required=True, type=Path)
    parser.add_argument("--leafcutter-counts", required=True, type=Path)
    parser.add_argument("--majiq-raw", required=True, type=Path)
    parser.add_argument("--tilgner-matrix", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--tealeaf-replication", type=Path)
    parser.add_argument("--paired-clr-input", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix, features, columns = read_tilgner_matrix(args.tilgner_matrix, args.gtf)
    mapped = set(columns["tealeaf_cell_type"].dropna())
    all_comparators = load_comparators(args.comparator_tests, mapped)
    selected = load_comparators(args.comparator_tests, mapped, significant_only=True)
    bundle = JunctionBundle.load(args.junction_bundle)
    groups = {}
    groups.update(scquint_groups(selected, bundle))
    groups.update(leafcutter_groups(selected, bundle, args.leafcutter_counts))
    groups.update(majiq_groups(selected, args.majiq_raw, bundle))
    if args.paired_clr_input:
        paired = pd.read_csv(args.paired_clr_input, sep="\t")
        paired_groups = scquint_groups(paired.assign(method="scQuint"), bundle)
    else:
        paired, paired_groups = add_paired_clr_discoveries(all_comparators, bundle, jobs=args.jobs)
    print(f"Selected {len(paired_groups):,} paired CLR discoveries", file=sys.stderr, flush=True)
    for group in paired_groups.values():
        group["method"] = "Paired junction CLR"
    groups.update({("Paired junction CLR", key[1], key[2]): value for key, value in paired_groups.items()})
    groups = add_short_read_effects(groups, bundle)
    coordinates = [(str(row.chromosome), int(row.start), int(row.end), str(row.strand)) for row in bundle.junctions.itertuples(index=False)]
    for group in groups.values():
        group["coordinates"] = [coordinates[index] for index in group["indices"]]
    results = assess_groups(groups, source_junction_counts(groups, matrix, features, columns, args.gtf))
    results.to_csv(args.output_dir / "junction_replication.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary = pd.concat([summarize(results), summarize(results, one_per_gene_contrast=True)], ignore_index=True)
    summary.to_csv(args.output_dir / "junction_replication_summary.tsv", sep="\t", index=False, na_rep="NA")
    if args.tealeaf_replication:
        combined = pd.concat([results, load_tealeaf_results(args.tealeaf_replication)], ignore_index=True, sort=False)
        comparison = pd.concat([summarize(combined), summarize(combined, one_per_gene_contrast=True)], ignore_index=True)
        comparison.to_csv(args.output_dir / "method_comparison_summary.tsv", sep="\t", index=False, na_rep="NA")
    paired.to_csv(args.output_dir / "paired_clr_discoveries.tsv.gz", sep="\t", index=False, na_rep="NA")
    manifest = {"source": "Joglekar et al. 2024 ScISOr-Seq2 processed annotated-transcript UMI matrix", "source_doi": "10.1038/s41593-024-01616-4", "selection": "full-data BH FDR below 0.05 within each method and cell-type contrast", "methods": ["LeafCutter", "MAJIQ Heterogen", "scQuint", "Paired junction CLR"], "junction_effect_estimator": "difference between equal-pseudobulk mean short-read junction compositions", "source_junction_estimator": "source transcript UMIs summed over annotated transcripts containing each junction", "eligibility": "every tested junction maps to at least one source annotated transcript; minimum pooled depth 20 per cell type; strict minimum depth 10 per cell type and biological replicate", "replication": "positive dot product between short-read and source long-read junction-composition differences", "paired_clr_minimum_subjects": 8, "paired_clr_bh_scope": "within cell-type contrast"}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
