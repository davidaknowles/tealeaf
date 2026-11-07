"""Summarize the model comparison: split replication, long-read direction, discoveries and primer offsets."""
import gzip, json
import numpy as np, pandas as pd
from tealeaf.sc.junction_benchmark import benjamini_hochberg
from tealeaf.sc.local_path_reads import block_structure, structure_paths

blocks = json.load(gzip.open(".cache/native_v4/collection_parity_p2_v4/all/blocks.json.gz", "rt"))
internal = {b["block_id"] for b in blocks.values() if block_structure(structure_paths(b))["structure"] == "internal"}
rows = []
for V in ("v4", "v4sh", "v4ip"):
    d = f".cache/internal_priming/compare/{V}"
    m = pd.read_csv(f"{d}/split_metrics.tsv", sep="\t")
    m = m.loc[~m.comparison.str.contains("hybrid")]
    tealeaf = m.loc[m.method.eq("Tealeaf")]
    out = {"model": V}
    for r in tealeaf.itertuples():
        out[f"split {r.comparison} {r.universe}"] = r.replicated_bh
    lr = pd.read_csv(f"{d}/lr_rank_summary.tsv", sep="\t")
    lr = lr.loc[lr.scope.eq("all tested") & lr.ranking.eq("continuous raw")]
    for r in lr.itertuples():
        out[f"LR A{r.cutoff}"] = round(r.normalized_auc, 3); out[f"LR agree{r.cutoff}"] = r.n_agree
    out["LR evaluable"] = int(lr.n_evaluable.iloc[0]) if "n_evaluable" in lr else np.nan
    fits = []
    for model, structure in ((f"internal_native_{V}", "internal"), (f"tj_native_{V}", "terminal")):
        f = pd.read_csv(f".cache/binary_grid/{model}_full/merged/paired_path.tsv", sep="\t", usecols=["test_id", "block_id", "gene_id", "level_a", "level_b", "p_value", "converged", "primer_offsets", "n_paths"])
        keep = f.block_id.isin(internal) if structure == "internal" else ~f.block_id.isin(internal)
        fits.append(f.loc[keep].assign(structure=structure))
    f = pd.concat(fits, ignore_index=True)
    p = f.p_value.where(f.converged.astype(bool), 1.).fillna(1.)
    f["fdr"] = benjamini_hochberg(p.to_numpy())
    disc = f.loc[f.fdr.le(.05)]
    out["full tests"] = len(f); out["converged"] = int(f.converged.astype(bool).sum())
    out["discoveries"] = len(disc); out["gene-pair discoveries"] = disc.drop_duplicates(["gene_id", "level_a", "level_b"]).shape[0]
    out["discoveries internal"] = int(disc.structure.eq("internal").sum()); out["discoveries terminal"] = int(disc.structure.eq("terminal").sum())
    offsets = np.concatenate([np.asarray(json.loads(x), float) for x in f.primer_offsets.dropna() if x.startswith("[")])
    out["offsets at |4|"] = round(float((np.abs(offsets) >= 3.99).mean()), 3); out["offset |median|"] = round(float(np.median(np.abs(offsets))), 3)
    rows.append(out)
t = pd.DataFrame(rows).set_index("model").T
t.to_csv(".cache/internal_priming/compare/summary.tsv", sep="\t")
print(t.to_string())
