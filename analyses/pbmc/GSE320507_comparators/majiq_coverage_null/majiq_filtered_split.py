"""Split replication with MAJIQ tests filtered by the coverage-only null (majiq_coverage_null.py on each fold).

Variants of the MAJIQ rows of a comparator split universe (build_split_universe.py):
  all             unchanged (reproduces the reported counts)
  null_p_ge_05    drop tests whose coverage-only null p (median over simulations) is < .05
  cov_ratio_le_05 keep tests with |log2 median LSV coverage ratio| <= .5 between the levels
  cov_ratio_le_1  keep tests with |log2 ratio| <= 1
Keys (gene, pair) must keep a test in both folds. Tealeaf fold fits are scored on the same universes.
Usage: majiq_filtered_split.py pbmc <GSE> <universe> <null dir prefix> <contrasts prefix> <output>
       majiq_filtered_split.py brain - <universe> <null dir prefix> <contrasts prefix> <output>
(<null dir prefix>{0,1}/majiq_coverage_null.parquet, <contrasts prefix>{0,1}/contrasts.json)
"""
import gzip, json, sys
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from extra_scripts.assess_paired_inference_audit import fit_table, split_universe_assessment  # noqa: E402
from tealeaf.sc.local_path_reads import block_structure, structure_paths  # noqa: E402

mode, gse, universe_path, null_prefix, contrasts_prefix, out = sys.argv[1:7]
out = Path(out)
out.mkdir(parents=True, exist_ok=True)
cache = REPO / ".cache/binary_grid"
if mode == "pbmc":
    folds = [pd.concat([fit_table(cache / f"pbmc_{gse}_{s}_{k}/merged/paired_path.tsv").assign(structure=s) for s in ("internal", "terminal")], ignore_index=True) for k in (0, 1)]
    name = "Tealeaf"
else:
    blocks = json.load(gzip.open(REPO / ".cache/native_v4/collection_parity_p2_v4/all/blocks.json.gz", "rt"))
    internal = {b["block_id"] for b in blocks.values() if block_structure(structure_paths(b))["structure"] == "internal"}
    folds = []
    for k in (0, 1):
        a = fit_table(cache / f"internal_native_v4_{k}/merged/paired_path.tsv")
        t = fit_table(cache / f"tj_native_v4_{k}/merged/paired_path.tsv")
        folds.append(pd.concat([a.loc[a.block_id.isin(internal)], t.loc[~t.block_id.isin(internal)]], ignore_index=True))
    name = "Tealeaf"

universe = pd.read_csv(universe_path, sep="\t", dtype={"gene_id": str}, low_memory=False)
majiq = universe.loc[universe.comparison.eq("MAJIQ Heterogen")].copy()
parts = []
for k in (0, 1):
    contrasts = json.loads(Path(f"{contrasts_prefix}{k}/contrasts.json").read_text())
    pair = {c["contrast_id"]: "||".join(sorted([c["level_a"], c["level_b"]])) for c in contrasts}
    null = pd.read_parquet(f"{null_prefix}{k}/majiq_coverage_null.parquet")
    null = null.dropna(subset=["feature_id"])
    null["pair_id"] = null.contrast_id.map(pair)
    null["p_null"] = null[[c for c in null if c.startswith("p_null")]].median(1)
    null = null[["pair_id", "feature_id", "p_null", "log2_cov_ratio"]].drop_duplicates(["pair_id", "feature_id"])
    local = majiq.loc[majiq.fold.eq(k)].merge(null, on=["pair_id", "feature_id"], how="left")
    print(f"fold {k}: MAJIQ rows {len(local)}, joined to null {int(local.p_null.notna().sum())}", flush=True)
    parts.append(local)
majiq = pd.concat(parts, ignore_index=True)
variants = {"all": np.ones(len(majiq), bool),
            "null_p_ge_05": ~(majiq.p_null < .05).to_numpy(),
            "cov_ratio_le_05": (majiq.log2_cov_ratio.abs() <= .5).to_numpy(),
            "cov_ratio_le_1": (majiq.log2_cov_ratio.abs() <= 1).to_numpy()}
frames, counts = [], []
for variant, mask in variants.items():
    local = majiq.loc[mask].copy()
    both = local.groupby(["gene_id", "pair_id"]).fold.nunique().pipe(lambda s: s[s == 2]).index
    local = local.set_index(["gene_id", "pair_id"]).loc[lambda d: d.index.isin(both)].reset_index()
    local["comparison"] = f"MAJIQ Heterogen [{variant}]"
    local["method"] = local["comparison"]
    counts.append(dict(variant=variant, tests=len(local), keys=len(both), genes=local.gene_id.nunique()))
    frames.append(local[universe.columns])
filtered = pd.concat(frames, ignore_index=True)
filtered.to_csv(out / "majiq_filtered_universe.tsv.gz", sep="\t", index=False)
pd.DataFrame(counts).to_csv(out / "universe_sizes.tsv", sep="\t", index=False)
print(pd.DataFrame(counts).to_string(index=False), flush=True)
split_universe_assessment(folds, out / "majiq_filtered_universe.tsv.gz", out, name)
