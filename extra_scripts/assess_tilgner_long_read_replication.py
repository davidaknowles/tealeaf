#!/usr/bin/env python3
"""Assess Tealeaf cell-type effects in Tilgner mouse-brain long reads."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import binomtest, chi2_contingency


ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]+)"')


def stable_identifier(value):
    return str(value).split(".", 1)[0]


def read_transcript_name_map(gtf_path):
    """Map stable gene ID and transcript name to a GENCODE transcript ID."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    mapping = {}
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            if not {"gene_id", "transcript_id", "transcript_name"} <= attributes.keys():
                continue
            key = (stable_identifier(attributes["gene_id"]), attributes["transcript_name"])
            mapping[key] = stable_identifier(attributes["transcript_id"])
    return mapping


def tilgner_cell_type(age, region, cell_type, subtype, oligodendrocyte_scope="mature"):
    """Map a Tilgner P56 annotation to the corresponding Tealeaf cell type."""
    if age != "P56":
        return None
    if cell_type == "Astro" and subtype in {"Astrocytes", "BergmannGlia"}:
        return "ASC"
    oligo_subtypes = {"MFOLs", "MOLs"} if oligodendrocyte_scope == "mature" else {"COPs", "DivOPCs", "MFOLs", "MOLs", "OPCs"}
    if cell_type == "Oligo" and subtype in oligo_subtypes:
        return "ODC"
    if region == "Cerebellum" and cell_type == "ExciteNeuron" and subtype == "GranuleCells":
        return "EX_cerebellum_granule"
    if region == "VisCortex" and cell_type == "ExciteNeuron" and subtype in {"ExciteL23", "ExciteL4", "ExciteL5", "ExciteL6"}:
        return "EX_cortical"
    if region == "Hippocampus" and cell_type == "ExciteNeuron" and subtype == "ExciteDG":
        return "EX_hippocampus_granule"
    if region == "Hippocampus" and cell_type == "ExciteNeuron" and subtype == "ExciteCA":
        return "EX_hippocampus_pyramidal"
    if region == "Thalamus" and cell_type == "ExciteNeuron" and subtype == "ExciteNeuron":
        return "EX_thalamus"
    if region == "Striatum" and cell_type == "InhibNeuron" and subtype in {"D1MSN", "D2MSN"}:
        return "INH_medium_spiny"
    return None


def biological_replicate(region, sample):
    """Return the replicate number used by the source study for P56 samples."""
    pairs = {
        "Cerebellum": {"M1": 1, "M2": 2},
        "Hippocampus": {"M1": 1, "M2": 2},
        "Striatum": {"M1": 1, "M2": 2},
        "Thalamus": {"M5": 1, "M6": 2},
        "VisCortex": {"M8": 1, "M9": 2},
    }
    return pairs.get(region, {}).get(sample)


def read_tilgner_matrix(matrix_dir, gtf_path, oligodendrocyte_scope="mature"):
    """Read annotated-transcript UMI counts and mapped P56 cell-type columns."""
    matrix_dir = Path(matrix_dir)
    features = pd.read_csv(matrix_dir / "features.tsv.gz", sep="\t", header=None, names=["gene_id", "transcript_name"])
    barcodes = pd.read_csv(matrix_dir / "barcodes.tsv.gz", sep="\t", header=None, names=["label"])
    parsed = barcodes["label"].str.split("::", expand=True)
    parsed.columns = ["age", "sample", "region", "broad", "cell_type", "subtype"]
    parsed["tealeaf_cell_type"] = [tilgner_cell_type(row.age, row.region, row.cell_type, row.subtype, oligodendrocyte_scope) for row in parsed.itertuples(index=False)]
    parsed["replicate"] = [biological_replicate(row.region, row.sample) for row in parsed.itertuples(index=False)]
    parsed["column"] = np.arange(len(parsed))
    transcript_names = read_transcript_name_map(gtf_path)
    features["stable_gene_id"] = features["gene_id"].map(stable_identifier)
    features["transcript_id"] = [transcript_names.get((gene, name)) for gene, name in zip(features["stable_gene_id"], features["transcript_name"])]
    features["row"] = np.arange(len(features))
    matrix = mmread(matrix_dir / "matrix.mtx.gz", spmatrix=True).tocsr()
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(f"matrix shape {matrix.shape} does not match features and barcodes")
    return matrix, features, parsed


def load_blocks(path):
    with gzip.open(path, "rt") as handle:
        blocks = json.load(handle)
    return {row["block_id"]: row for row in blocks}


def load_path_usage(root, selected_tests):
    tables = []
    for path in sorted(Path(root).glob("shard_*/path_usage.tsv")):
        table = pd.read_csv(path, sep="\t", usecols=["test_id", "cell_type", "path_number", "proportion"])
        table = table[table["test_id"].isin(selected_tests)]
        if not table.empty:
            tables.append(table)
    if not tables:
        raise FileNotFoundError(f"no selected path usages under {root}")
    return pd.concat(tables, ignore_index=True).groupby(["test_id", "cell_type", "path_number"], as_index=False)["proportion"].mean()


def block_feature_rows(block, tested_signatures, features):
    signature_to_path = {tuple(tuple(interval) for interval in signature): index + 1 for index, signature in enumerate(tested_signatures)}
    transcript_paths = {}
    for transcript, path_index in zip(block["transcripts"], block["path_index"]):
        signature = tuple(tuple(interval) for interval in block["path_signatures"][path_index])
        if signature in signature_to_path:
            transcript_paths[stable_identifier(transcript)] = signature_to_path[signature]
    gene = stable_identifier(block["gene_id"])
    local = features[(features["stable_gene_id"] == gene) & features["transcript_id"].notna()].copy()
    local["path_number"] = local["transcript_id"].map(transcript_paths)
    return local[local["path_number"].notna()].astype({"path_number": int})


def normalized_difference(first, second):
    first_total, second_total = first.sum(), second.sum()
    if first_total <= 0 or second_total <= 0:
        return np.full(len(first), np.nan)
    return second / second_total - first / first_total


def vector_agreement(original, external):
    original_norm = np.linalg.norm(original)
    external_norm = np.linalg.norm(external)
    if original_norm <= 0 or external_norm <= 0:
        return np.nan, np.nan
    dot = float(original @ external)
    return dot, dot / (original_norm * external_norm)


def bh_adjust(p_values):
    values = np.asarray(p_values, dtype=float)
    output = np.full(len(values), np.nan)
    valid = np.isfinite(values)
    order = np.argsort(values[valid])
    ranked = values[valid][order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    valid_positions = np.flatnonzero(valid)[order]
    output[valid_positions] = np.minimum(adjusted, 1.0)
    return output


def wilson_interval(successes, total, z=1.959963984540054):
    if total == 0:
        return np.nan, np.nan
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half_width = z * np.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return center - half_width, center + half_width


def summarize_replication(results):
    rows = []
    for minimum_depth in (10, 20, 50, 100):
        for minimum_effect in (0.0, 0.05, 0.1, 0.2, 0.3):
            pooled = results[(results["mapping_complete"]) & (results["minimum_pooled_depth"] >= minimum_depth) & (results["original_effect_norm"] >= minimum_effect) & results["pooled_replicated"].notna()]
            strict = pooled[pooled["minimum_replicate_depth"] >= minimum_depth / 2]
            for endpoint, group, column in (("pooled direction", pooled, "pooled_replicated"), ("both biological replicates", strict, "both_replicates_replicated")):
                values = group[column].dropna().astype(bool)
                successes = int(values.sum())
                low, high = wilson_interval(successes, len(values))
                local_fdr = bh_adjust(group.loc[values.index, "technical_g_test_p_value"])
                technical_directional = (local_fdr < 0.05) & values.to_numpy()
                if endpoint == "pooled direction":
                    null_rate = 0.5
                    orientation_trials = len(values)
                else:
                    replicate_signs_agree = ((group.loc[values.index, "replicate_1_dot_product"] > 0) == (group.loc[values.index, "replicate_2_dot_product"] > 0)).to_numpy()
                    orientation_trials = int(replicate_signs_agree.sum())
                    null_rate = orientation_trials / (2 * len(values)) if len(values) else np.nan
                rows.append({"endpoint": endpoint, "minimum_depth": minimum_depth, "minimum_original_effect_norm": minimum_effect, "n_tests": len(values), "n_replicated": successes, "replication_rate": successes / len(values) if len(values) else np.nan, "ci_low": low, "ci_high": high, "null_rate": null_rate, "n_orientation_trials": orientation_trials, "sign_test_p_value": binomtest(successes, orientation_trials, 0.5, alternative="greater").pvalue if orientation_trials else np.nan, "n_technical_bh_directional": int(technical_directional.sum()), "technical_bh_directional_rate": technical_directional.mean() if len(values) else np.nan})
    return pd.DataFrame(rows)


def assess_replication(catalog, usage, blocks, matrix, features, columns):
    column_groups = {}
    for (cell_type, replicate), group in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby(["tealeaf_cell_type", "replicate"]):
        column_groups[(cell_type, int(replicate))] = group["column"].to_numpy(dtype=int)
    usage_lookup = usage.set_index(["test_id", "cell_type", "path_number"])["proportion"]
    rows = []
    for record in catalog.itertuples(index=False):
        parts = record.test_id.split("|")
        level_a, level_b = parts[-2:]
        signatures = json.loads(record.tested_path_signatures)
        n_paths = len(signatures)
        block = blocks.get(record.block_id)
        local_features = block_feature_rows(block, signatures, features) if block is not None else pd.DataFrame()
        row_groups = {path: group["row"].to_numpy(dtype=int) for path, group in local_features.groupby("path_number")}
        original_a = np.array([usage_lookup.get((record.test_id, level_a, path), np.nan) for path in range(1, n_paths + 1)])
        original_b = np.array([usage_lookup.get((record.test_id, level_b, path), np.nan) for path in range(1, n_paths + 1)])
        original_delta = original_b - original_a
        counts = {}
        for level in (level_a, level_b):
            for replicate in (1, 2):
                selected_columns = column_groups.get((level, replicate), np.array([], dtype=int))
                counts[(level, replicate)] = np.array([matrix[rows_for_path][:, selected_columns].sum() if len(selected_columns) and len(rows_for_path) else 0.0 for rows_for_path in (row_groups.get(path, np.array([], dtype=int)) for path in range(1, n_paths + 1))], dtype=float)
        pooled_a = counts[(level_a, 1)] + counts[(level_a, 2)]
        pooled_b = counts[(level_b, 1)] + counts[(level_b, 2)]
        pooled_delta = normalized_difference(pooled_a, pooled_b)
        pooled_dot, pooled_cosine = vector_agreement(original_delta, pooled_delta)
        replicate_dots, replicate_cosines = [], []
        for replicate in (1, 2):
            delta = normalized_difference(counts[(level_a, replicate)], counts[(level_b, replicate)])
            dot, cosine = vector_agreement(original_delta, delta)
            replicate_dots.append(dot)
            replicate_cosines.append(cosine)
        contingency = np.vstack([pooled_a, pooled_b])
        nonzero = contingency.sum(axis=0) > 0
        technical_p = np.nan
        if nonzero.sum() >= 2 and np.all(contingency[:, nonzero].sum(axis=1) > 0):
            technical_p = chi2_contingency(contingency[:, nonzero], lambda_="log-likelihood").pvalue
        mapping_complete = len(row_groups) == n_paths and all((level, replicate) in column_groups for level in (level_a, level_b) for replicate in (1, 2)) and np.isfinite(original_delta).all()
        rows.append({"test_id": record.test_id, "block_id": record.block_id, "gene_id": record.gene_id, "gene_name": record.gene_name, "event_type": record.event_type, "level_a": level_a, "level_b": level_b, "n_paths": n_paths, "n_mappable_paths": len(row_groups), "n_mapped_transcripts": len(local_features), "mapping_complete": mapping_complete, "original_effect_norm": float(np.linalg.norm(original_delta)) if np.isfinite(original_delta).all() else np.nan, "long_read_effect_norm": float(np.linalg.norm(pooled_delta)) if np.isfinite(pooled_delta).all() else np.nan, "pooled_dot_product": pooled_dot, "pooled_cosine": pooled_cosine, "pooled_replicated": pooled_dot > 0 if np.isfinite(pooled_dot) else np.nan, "replicate_1_dot_product": replicate_dots[0], "replicate_2_dot_product": replicate_dots[1], "replicate_1_cosine": replicate_cosines[0], "replicate_2_cosine": replicate_cosines[1], "both_replicates_replicated": replicate_dots[0] > 0 and replicate_dots[1] > 0 if np.isfinite(replicate_dots).all() else np.nan, "minimum_pooled_depth": float(min(pooled_a.sum(), pooled_b.sum())), "minimum_replicate_depth": float(min(*(counts[(level, replicate)].sum() for level in (level_a, level_b) for replicate in (1, 2)))), "technical_g_test_p_value": technical_p, "counts_a_rep1": json.dumps(counts[(level_a, 1)].astype(int).tolist()), "counts_a_rep2": json.dumps(counts[(level_a, 2)].astype(int).tolist()), "counts_b_rep1": json.dumps(counts[(level_b, 1)].astype(int).tolist()), "counts_b_rep2": json.dumps(counts[(level_b, 2)].astype(int).tolist())})
    results = pd.DataFrame(rows)
    results["technical_g_test_fdr"] = bh_adjust(results["technical_g_test_p_value"])
    results["technical_bh_directional_replication"] = (results["technical_g_test_fdr"] < 0.05) & (results["pooled_replicated"] == True)
    return results


def write_plot(summary, output_path):
    from plotnine import aes, coord_cartesian, element_blank, element_text, facet_wrap, geom_col, geom_errorbar, geom_hline, geom_text, ggplot, labs, position_dodge, scale_fill_manual, theme, theme_bw
    selected = summary[(summary["minimum_depth"] == 20) & (summary["minimum_original_effect_norm"].isin([0.0, 0.05, 0.1, 0.2]))].copy()
    selected["effect threshold"] = selected["minimum_original_effect_norm"].map(lambda value: f"at least {value:g}")
    selected["label"] = selected.apply(lambda row: f"{row.replication_rate:.0%}\n({int(row.n_tests)})" if row.n_tests else "NA", axis=1)
    baselines = selected[["endpoint", "null_rate"]].drop_duplicates()
    plot = ggplot(selected, aes("effect threshold", "replication_rate", fill="endpoint"))
    plot += geom_hline(data=baselines, mapping=aes(yintercept="null_rate"), inherit_aes=False, linetype="dashed", color="#666666")
    plot += geom_col(position=position_dodge(width=0.8), width=0.72)
    plot += geom_errorbar(aes(ymin="ci_low", ymax="ci_high"), position=position_dodge(width=0.8), width=0.16)
    plot += geom_text(aes(label="label"), position=position_dodge(width=0.8), va="bottom", size=8, nudge_y=0.02)
    plot += facet_wrap("~ endpoint")
    plot += scale_fill_manual(values={"pooled direction": "#0B6666", "both biological replicates": "#A94E16"})
    plot += coord_cartesian(ylim=(0, 1))
    plot += labs(x="Minimum Tealeaf path-effect norm", y="Directional replication rate", title="Independent long reads support Tealeaf cell-type path effects", caption="Dashed lines show the conditional sign-flip null:\n50% pooled; half the observed between-replicate concordance for the strict endpoint.")
    plot += theme_bw(base_size=10)
    plot += theme(legend_position="none", panel_grid_minor=element_blank(), axis_text_x=element_text(rotation=30, ha="right"), plot_title=element_text(size=11), plot_caption=element_text(size=8))
    plot.save(output_path, width=8.3, height=4.8, units="in", verbose=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tilgner-matrix", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--block-cache", required=True, type=Path)
    parser.add_argument("--pairwise-catalog", required=True, type=Path)
    parser.add_argument("--path-usage-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--oligodendrocyte-scope", choices=["mature", "all"], default="mature")
    parser.add_argument("--source-commit")
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(args.pairwise_catalog, sep="\t")
    catalog = catalog[pd.to_numeric(catalog["fdr"], errors="coerce") < 0.05].copy()
    usage = load_path_usage(args.path_usage_root, set(catalog["test_id"]))
    matrix, features, columns = read_tilgner_matrix(args.tilgner_matrix, args.gtf, args.oligodendrocyte_scope)
    results = assess_replication(catalog, usage, load_blocks(args.block_cache), matrix, features, columns)
    summary = summarize_replication(results)
    results.to_csv(args.output_dir / "pairwise_replication.tsv", sep="\t", index=False, na_rep="NA")
    summary.to_csv(args.output_dir / "replication_summary.tsv", sep="\t", index=False, na_rep="NA")
    columns.to_csv(args.output_dir / "tilgner_column_mapping.tsv", sep="\t", index=False, na_rep="NA")
    manifest = {"source": "Joglekar et al. 2024 ScISOr-Seq2 processed annotated-transcript UMI matrix", "source_doi": "10.1038/s41593-024-01616-4", "source_repository": "git@github.com:noush-joglekar/biccn_tilgner_scisorseq.git", "source_commit": args.source_commit, "age": "P56", "oligodendrocyte_scope": args.oligodendrocyte_scope, "statistical_selection": "Tealeaf pairwise BH FDR below 0.05 at A_test=64", "effect_estimates": "Tealeaf A_report=1 subject-mean path proportions", "replication_definition": "positive dot product between Tealeaf and long-read path-proportion difference vectors", "strict_replication_definition": "positive dot product independently in both Tilgner biological replicates", "strict_null_definition": "half the observed between-replicate directional concordance under randomized Tealeaf contrast orientation", "cell_type_mapping": "fixed mapping implemented in extra_scripts/assess_tilgner_long_read_replication.py", "technical_test": "pooled long-read UMI G-test, provided as a secondary read-level endpoint rather than biological replication"}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not args.skip_plot:
        write_plot(summary, args.output_dir / "long_read_replication.pdf")


if __name__ == "__main__":
    main()
