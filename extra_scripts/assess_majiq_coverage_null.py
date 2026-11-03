#!/usr/bin/env python3
"""Coverage-only null of MAJIQ Heterogen tests (tealeaf.sc.majiq_coverage) for a set of contrasts.

Reads per-sample LSV totals, connection reads and quantifiability from the psicov store of a MAJIQ build, runs
coverage_null_contrast for every contrast with a heterogen raw table, and joins each test to that table (normalized
feature_id, MAJIQ's raw p-value). Writes majiq_coverage_null.parquet and summaries: validation.tsv (agreement of the
reproduced p with MAJIQ's), null_calls.tsv, by_pooled_effect.tsv, by_coverage_ratio.tsv. Needs the rna_majiq package.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.junction_benchmark import contrast_file_stem
from tealeaf.sc.majiq_coverage import coverage_null_contrast

KEY = ["gene_id", "event_type", "ref_exon_start", "ref_exon_end", "start", "end"]
FEATURE = ["gene_id", "seqid", "strand", "event_type", "ref_exon_start", "ref_exon_end", "start", "end", "is_intron", "other_exon_start", "other_exon_end"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True, help="directory with all.psicov and combined.splicegraph")
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True, help="heterogen raw tables <contrast stem>.tsv")
    parser.add_argument("--min-experiments", type=float, required=True, help="per group, as given to heterogen; < 1 is a fraction of the group")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n-sim", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    import rna_majiq as nm

    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    psicov = nm.PsiCoverage.from_zarr(args.build / "all.psicov")
    prefixes = {sample: index for index, sample in enumerate(psicov.prefixes)}
    total = np.asarray(psicov.df["raw_total"].values, dtype=float)
    cover = np.rint(np.nan_to_num(np.asarray(psicov.df["raw_psi"].values, dtype=float)) * total)
    passed = np.asarray(psicov.df["event_passed"].values)
    alpha = np.asarray(psicov.alpha_prior.values, dtype=float)
    graph = nm.SpliceGraph.from_zarr(args.build / "combined.splicegraph")
    connections = psicov.get_events(graph.introns, graph.junctions).ec_dataframe().reset_index()[["ec_idx"] + KEY]
    connections["event_type"] = connections.event_type.map(lambda value: value.decode() if isinstance(value, bytes) else str(value))
    tables = []
    for contrast in json.loads(args.contrasts.read_text()):
        raw_path = args.raw_dir / f"{contrast_file_stem(contrast['contrast_id'])}.tsv"
        if not raw_path.exists():
            continue
        columns = np.array([prefixes[s] for s in contrast["samples_a"] if s in prefixes] + [prefixes[s] for s in contrast["samples_b"] if s in prefixes])
        in_a = np.arange(len(columns)) < sum(s in prefixes for s in contrast["samples_a"])
        table = coverage_null_contrast(total[:, columns], cover[:, columns], passed[:, columns], alpha, in_a, args.min_experiments, rng, args.n_sim)
        if table.empty:
            continue
        table.insert(0, "contrast_id", contrast["contrast_id"])
        table = table.merge(connections, on="ec_idx", how="left")
        raw = pd.read_csv(raw_path, sep="\t", comment="#", low_memory=False)
        raw["event_type"] = raw.event_type.astype(str)
        raw["feature_id"] = raw[FEATURE].astype(str).agg(":".join, axis=1)
        tables.append(table.merge(raw[KEY + ["feature_id", "raw_pvalue"]], on=KEY, how="left"))
        print(contrast["contrast_id"], len(table), flush=True)
    result = pd.concat(tables, ignore_index=True)
    result.drop(columns=KEY).to_parquet(args.output / "majiq_coverage_null.parquet")
    summarize(result, args.n_sim, args.output)


def summarize(table, n_sim, output):
    joined = table.dropna(subset=["raw_pvalue"])
    validation = {"tests": len(table), "joined": len(joined), "spearman_p": float(joined[["p_obs", "raw_pvalue"]].corr("spearman").iloc[0, 1]),
                  "median_abs_log10_ratio": float(np.nanmedian(np.abs(np.log10(joined.p_obs / joined.raw_pvalue))))}
    call = table.q_obs < .05
    calls = pd.DataFrame([{"sim": k, "tests": len(table), "observed_q05": int(call.sum()), "null_p05": float((table[f"p_null{k}"] < .05).mean()),
                           "null_q05": int((table[f"q_null{k}"] < .05).sum())} for k in range(n_sim)])
    null_call = table.q_null0 < .05
    same = np.sign(table.majiq_effect) == np.sign(table[[f"effect_null{k}" for k in range(n_sim)]].mean(axis=1))
    by_effect = []
    for label, group in table.groupby(pd.cut(table.pooled_effect.abs(), [-1e-9, .02, .05, .1, .2, 1]), observed=True):
        mask = call.loc[group.index]
        by_effect.append({"pooled_effect": str(label), "tests": len(group), "observed_q05": int(mask.sum()), "null_q05": int(null_call.loc[group.index].sum()),
                          "null_p05_rate": float((group.p_null0 < .05).mean()), "calls_also_null_q05": int((mask & null_call.loc[group.index]).sum()),
                          "calls_sign_as_null": float(same.loc[group.index][mask].mean())})
    by_ratio = [{"abs_log2_coverage_ratio": str(label), "tests": len(group), "observed_q05_rate": float(call.loc[group.index].mean()),
                 "null_q05_rate": float(null_call.loc[group.index].mean())}
                for label, group in table.groupby(pd.cut(table.log2_cov_ratio.abs(), [-1e-9, .5, 1, 2, 3, 100]), observed=True)]
    pd.DataFrame([validation]).to_csv(output / "validation.tsv", sep="\t", index=False)
    calls.to_csv(output / "null_calls.tsv", sep="\t", index=False)
    pd.DataFrame(by_effect).to_csv(output / "by_pooled_effect.tsv", sep="\t", index=False)
    pd.DataFrame(by_ratio).to_csv(output / "by_coverage_ratio.tsv", sep="\t", index=False)
    print(validation, calls.to_string(index=False), pd.DataFrame(by_effect).to_string(index=False), pd.DataFrame(by_ratio).to_string(index=False), sep="\n", flush=True)


if __name__ == "__main__":
    main()
