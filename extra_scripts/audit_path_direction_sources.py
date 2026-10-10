#!/usr/bin/env python3
"""Compare EC and junction-only local path directions with long reads.

For every full-data pairwise local-path test, the EC direction is the existing
weak-prior (zeta=1) subject-mean path-proportion difference. The junction
direction re-estimates the same block paths from STARsolo junction UMIs only.
Both are scored against the same long-read path effect, on the same tests and
under the same ranking, so a difference in agreement isolates the measurement.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, load_blocks, block_feature_rows, normalized_difference, vector_agreement
from tealeaf.sc.junction_benchmark import JunctionBundle
from tealeaf.sc.junction_paths import block_path_exons, path_junction_map, junction_identifiable, junction_path_usage
from tealeaf.sc.replication_audit import ranked_direction_summary


def read_shards(root, name, **kwargs):
    return pd.concat([pd.read_csv(path, sep="\t", **kwargs) for path in sorted(Path(root).glob(f"shard_*/{name}"))], ignore_index=True)


def anchor(value):
    value = json.loads(value) if isinstance(value, str) else value
    return None if value is None else tuple(value)


def junction_index(bundle):
    """(chromosome, STAR intron start, end) -> junction columns, all strands."""
    index = defaultdict(list)
    for column, row in enumerate(bundle.junctions[["chromosome", "start", "end"]].itertuples(index=False)):
        index[(row.chromosome, int(row.start), int(row.end))].append(column)
    return index


def junction_counts(bundle, index, chromosome, keys, sample_rows):
    """Subject x junction UMI matrix; block keys are (exon end, next exon start)."""
    columns = [index.get((chromosome, start + 1, end), []) for start, end in keys]
    found = np.array([len(item) > 0 for item in columns])
    counts = np.zeros((len(sample_rows), len(keys)))
    valid = [row for row in sample_rows if row is not None]
    if valid:
        local = bundle.counts[valid]
        for position, junction_columns in enumerate(columns):
            if junction_columns:
                values = np.asarray(local[:, junction_columns].sum(axis=1)).ravel()
                counts[[i for i, row in enumerate(sample_rows) if row is not None], position] = values
    return counts, found


def level_usage(counts, weights, concentration):
    pooled = counts.sum(axis=0)
    pooled_usage = junction_path_usage(pooled, weights, concentration) if pooled.sum() > 0 else np.full(weights.shape[1], np.nan)
    subject = [junction_path_usage(row, weights, concentration) for row in counts if row.sum() > 0]
    subject_usage = np.mean(subject, axis=0) if subject else np.full(weights.shape[1], np.nan)
    return pooled_usage, subject_usage, float(pooled.sum())


def assess(args):
    tests = read_shards(args.usage_root, "paired_path.tsv", usecols=["test_id", "block_id", "gene_id", "level_a", "level_b", "n_paths", "n_subjects", "converged", "median_gene_umis"])
    usage = read_shards(args.usage_root, "path_usage.tsv", usecols=["test_id", "subject", "cell_type", "path_number", "path_signature", "proportion"])
    ranking = pd.read_csv(args.ranking, sep="\t", usecols=["test_id", "p_value", "raw_p_value", "statistic", "fdr"])
    tests = tests.merge(ranking, on="test_id", how="left", validate="one_to_one")
    blocks = load_blocks(args.block_cache)
    bundle = JunctionBundle.load(args.junction_prefix)
    index = junction_index(bundle)
    sample_row = {(row.subject, row.cell_type): position for position, row in enumerate(bundle.samples.itertuples(index=False))}
    matrix, features, columns = read_tilgner_matrix(args.matrix_dir, args.gtf)
    groups = {(level, int(rep)): frame.column.to_numpy(int) for (level, rep), frame in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby(["tealeaf_cell_type", "replicate"])}
    feature_groups = {gene: frame for gene, frame in features.groupby("stable_gene_id")}
    usage_groups = dict(tuple(usage.groupby("test_id")))
    block_cache = {}
    records = []
    for row in tests.itertuples(index=False):
        local = usage_groups.get(row.test_id)
        block = blocks.get(row.block_id)
        if local is None or block is None:
            continue
        signatures = local.drop_duplicates("path_number").sort_values("path_number").path_signature.map(json.loads).tolist()
        n_paths = len(signatures)
        means = local.groupby(["cell_type", "path_number"]).proportion.mean()
        ec_a = np.array([means.get((row.level_a, path), np.nan) for path in range(1, n_paths + 1)])
        ec_b = np.array([means.get((row.level_b, path), np.nan) for path in range(1, n_paths + 1)])
        ec_effect = ec_b - ec_a
        key = (row.block_id, json.dumps(signatures))
        if key not in block_cache:
            paths = [block_path_exons(anchor(block["left_anchor"]), signature, anchor(block["right_anchor"])) for signature in signatures]
            keys, weights = path_junction_map(paths, block["strand"])
            gene = row.gene_id.split(".")[0]
            lr_rows = block_feature_rows(block, signatures, feature_groups.get(gene, features.iloc[:0]))
            path_rows = {path: frame.row.to_numpy(int) for path, frame in lr_rows.groupby("path_number")}
            complete = set(path_rows) == set(range(1, n_paths + 1))
            lr = {group: np.asarray([matrix[path_rows.get(path, np.empty(0, int))][:, cols].sum() for path in range(1, n_paths + 1)], dtype=float) for group, cols in groups.items()}
            block_cache[key] = keys, weights, complete, lr
        keys, weights, complete, lr = block_cache[key]
        lr_a = sum(lr.get((row.level_a, rep), np.zeros(n_paths)) for rep in (1, 2))
        lr_b = sum(lr.get((row.level_b, rep), np.zeros(n_paths)) for rep in (1, 2))
        lr_complete = complete and all((level, rep) in groups for level in (row.level_a, row.level_b) for rep in (1, 2))
        lr_effect = normalized_difference(lr_a, lr_b)
        lr_depth = min(lr_a.sum(), lr_b.sum())
        identifiable = junction_identifiable(weights)
        record = {"test_id": row.test_id, "block_id": row.block_id, "gene_id": row.gene_id, "level_a": row.level_a, "level_b": row.level_b, "n_paths": n_paths, "n_subjects": row.n_subjects, "median_gene_umis": row.median_gene_umis, "p_value": row.p_value, "raw_p_value": row.raw_p_value, "statistic": row.statistic, "fdr": row.fdr, "n_variable_junctions": len(keys), "junction_identifiable": identifiable, "lr_complete": lr_complete, "lr_depth": lr_depth, "lr_effect": json.dumps(lr_effect.tolist()), "lr_effect_norm": float(np.linalg.norm(lr_effect)) if np.isfinite(lr_effect).all() else np.nan, "ec_effect": json.dumps(ec_effect.tolist())}
        record["ec_dot"], _ = vector_agreement(ec_effect, lr_effect)
        junction = {"pooled": np.full(n_paths, np.nan), "subject": np.full(n_paths, np.nan)}
        depth, found = 0., 0
        if identifiable:
            subjects = sorted(set(local.subject))
            per_level = []
            for level in (row.level_a, row.level_b):
                counts, present = junction_counts(bundle, index, block["chromosome"], keys, [sample_row.get((subject, level)) for subject in subjects])
                per_level.append(level_usage(counts, weights, args.concentration))
                found = int(present.sum())
            junction = {"pooled": per_level[1][0] - per_level[0][0], "subject": per_level[1][1] - per_level[0][1]}
            depth = min(per_level[0][2], per_level[1][2])
        record.update({"junction_keys_found": found, "junction_depth": depth})
        for name, effect in junction.items():
            record[f"junction_{name}_effect"] = json.dumps(np.asarray(effect).tolist())
            record[f"junction_{name}_dot"], _ = vector_agreement(effect, lr_effect)
            record[f"ec_junction_{name}_dot"], _ = vector_agreement(ec_effect, effect)
        records.append(record)
    return pd.DataFrame(records)


def ec_effects(root):
    """test_id -> zeta=1 subject-mean path-proportion difference for one refit."""
    tests = read_shards(root, "paired_path.tsv", usecols=["test_id", "level_a", "level_b"])
    usage = read_shards(root, "path_usage.tsv", usecols=["test_id", "cell_type", "path_number", "proportion"])
    means = usage.groupby(["test_id", "cell_type", "path_number"]).proportion.mean()
    effects = {}
    present = set(means.index.get_level_values("test_id"))
    for row in tests.itertuples(index=False):
        if row.test_id not in present:
            continue
        local = means.loc[row.test_id]
        n_paths = int(local.index.get_level_values("path_number").max())
        values = [np.array([local.get((level, path), np.nan) for path in range(1, n_paths + 1)]) for level in (row.level_a, row.level_b)]
        effects[row.test_id] = values[1] - values[0]
    return effects


def add_ec_variant(table, root, name):
    """Score another EC refit against the stored junction and LR effects."""
    effects = ec_effects(root)
    dots, junction_dots = [], []
    for row in table.itertuples(index=False):
        effect = effects.get(row.test_id, np.array([np.nan]))
        lr_effect, junction = np.asarray(json.loads(row.lr_effect)), np.asarray(json.loads(row.junction_pooled_effect))
        dots.append(vector_agreement(effect, lr_effect)[0] if effect.shape == lr_effect.shape else np.nan)
        junction_dots.append(vector_agreement(effect, junction)[0] if effect.shape == junction.shape else np.nan)
    table[f"ec_{name}_dot"], table[f"ec_{name}_junction_dot"] = dots, junction_dots
    return table


def variant_summary(table, names):
    """Agreement of every EC variant and junctions on one common test set."""
    base = table.loc[table.lr_complete & table.lr_depth.ge(20) & np.isfinite(table.junction_pooled_dot) & table.junction_pooled_dot.ne(0)]
    rows = []
    for min_depth in (1, 20, 100):
        for min_lr in (0., .1):
            local = base.loc[base.junction_depth.ge(min_depth) & base.lr_effect_norm.ge(min_lr)]
            for name in names:
                local = local.loc[np.isfinite(local[f"ec_{name}_dot"]) & local[f"ec_{name}_dot"].ne(0)]
            row = {"min_junction_depth": min_depth, "min_lr_effect_norm": min_lr, "n_tests": len(local), "n_genes": local.gene_id.nunique(), "junction_pooled_agree_lr": local.junction_pooled_dot.gt(0).mean()}
            for name in names:
                row[f"ec_{name}_agree_lr"] = local[f"ec_{name}_dot"].gt(0).mean()
                row[f"ec_{name}_agree_junction"] = local[f"ec_{name}_junction_dot"].gt(0).mean()
            rows.append(row)
    return pd.DataFrame(rows)


def summarize(table):
    rows = []
    eligible = table.loc[table.lr_complete & table.lr_depth.ge(20) & np.isfinite(table.ec_dot) & table.ec_dot.ne(0)]
    for min_depth in (1, 20, 100, 500):
        for min_lr in (0., .1):
            local = eligible.loc[eligible.junction_depth.ge(min_depth) & eligible.lr_effect_norm.ge(min_lr)]
            for name in ("pooled", "subject"):
                usable = local.loc[np.isfinite(local[f"junction_{name}_dot"]) & local[f"junction_{name}_dot"].ne(0)]
                disagree = usable.loc[usable[f"ec_junction_{name}_dot"] < 0]
                rows.append({"min_junction_depth": min_depth, "min_lr_effect_norm": min_lr, "junction_estimator": name, "n_tests": len(usable), "n_genes": usable.gene_id.nunique(), "ec_agree_lr": usable.ec_dot.gt(0).mean(), "junction_agree_lr": usable[f"junction_{name}_dot"].gt(0).mean(), "ec_agree_junction": usable[f"ec_junction_{name}_dot"].gt(0).mean(), "n_ec_junction_disagree": len(disagree), "lr_sides_with_junction": disagree[f"junction_{name}_dot"].gt(0).mean() if len(disagree) else np.nan})
    return pd.DataFrame(rows)


def ranked(table):
    """Same tests and same production ordering, only the direction differs."""
    eligible = table.loc[table.lr_complete & table.lr_depth.ge(20) & np.isfinite(table.ec_dot) & table.ec_dot.ne(0) & np.isfinite(table.junction_pooled_dot) & table.junction_pooled_dot.ne(0)]
    rows = []
    for ranking, columns in (("calibrated then raw", ["p_value", "raw_p_value"]), ("continuous raw", ["raw_p_value", "p_value"])):
        ordered = eligible.sort_values(columns + ["statistic", "test_id"], ascending=[True, True, False, True], kind="stable")
        for name, column in (("EC subject mean", "ec_dot"), ("junction pooled", "junction_pooled_dot"), ("junction subject mean", "junction_subject_dot")):
            local = ordered.loc[np.isfinite(ordered[column]) & ordered[column].ne(0)].copy()
            local["rank"] = np.arange(1, len(local) + 1)
            local["method"] = name
            local["pooled_replicated"] = local[column].gt(0)
            rows.extend({**item, "ranking": ranking} for item in ranked_direction_summary(local))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--usage-root", type=Path, required=True)
    parser.add_argument("--ranking", type=Path, required=True)
    parser.add_argument("--junction-prefix", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--concentration", type=float, default=1.)
    parser.add_argument("--variant", action="append", default=[], help="name=refit root; scored against a stored table from --reuse-table")
    parser.add_argument("--reuse-table", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.reuse_table is not None:
        table = pd.read_csv(args.reuse_table, sep="\t")
        table["ec_production_dot"], table["ec_production_junction_dot"] = table.ec_dot, table.ec_junction_pooled_dot
        names = ["production"]
        for item in args.variant:
            name, root = item.split("=", 1)
            table = add_ec_variant(table, Path(root), name)
            names.append(name)
        summary = variant_summary(table, names)
        summary.to_csv(args.output_dir / "variant_summary.tsv", sep="\t", index=False, na_rep="NA")
        table.to_csv(args.output_dir / "variant_directions.tsv.gz", sep="\t", index=False, na_rep="NA")
        print(summary.to_string(index=False), flush=True)
        return
    table = assess(args)
    table.to_csv(args.output_dir / "path_direction_sources.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary = summarize(table)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    rank = ranked(table)
    rank.to_csv(args.output_dir / "rank_summary.tsv", sep="\t", index=False, na_rep="NA")
    print(f"tests {len(table)}, junction identifiable {table.junction_identifiable.sum()}, junction keys found {table.junction_keys_found.sum()} of {table.loc[table.junction_identifiable].n_variable_junctions.sum()}", flush=True)
    print(summary.to_string(index=False), flush=True)
    print(rank.to_string(index=False), flush=True)
    (args.output_dir / "manifest.json").write_text(json.dumps({"ec_direction": "zeta=1 subject-mean path-proportion difference from the existing full-data pairwise refit", "junction_direction": f"same block paths re-estimated from STARsolo junction UMIs only, EM with total concentration {args.concentration}", "long_read": "pooled Tilgner annotated-transcript UMIs collapsed by the same block path map, eligible at >=20 UMIs per cell type", "ranking": "production full-data zeta=64 pairwise p-values, identical for every direction estimator", "inputs": {key: str(value) for key, value in vars(args).items()}}, indent=2) + "\n")


if __name__ == "__main__":
    main()
