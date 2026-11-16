"""Junction-level view of rMATS's strongest plus-strand MXE calls in GSE320507 that have a nonzero long-read effect.

For each event (upstream exon U, exons 1 and 2, downstream exon D; fromGTF.MXE.txt), STARsolo junction counts summed over
the full-data pseudobulks of the two cell types: U-1 and 1-D (rMATS inclusion form), U-2 and 2-D (exclusion form), and 1-2
(both exons in one molecule, which counts toward neither form in the long reads but whose outer junctions U-1 and 2-D
count toward both rMATS forms). Reports the share of the event's junction reads on 1-2 and on the form-specific junctions
1-D and U-2, and the inclusion shares per cell type with all four form junctions and with 1-D and U-2 only.
Usage: inspect_mxe.py <output.tsv>"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
sys.path.insert(0, ".")
from extra_scripts.summarize_tilgner_comparator_ranks import evaluable, read_source_table

B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
lr = evaluable(read_source_table(Path("analyses/pbmc/GSE320507_comparators/lr_event/rmats_full_data_lr.tsv.gz")))
lr = lr.loc[lr.feature_id.str.startswith("rMATS:MXE:")].copy()
lr["id"] = lr.feature_id.str.split(":").str[2].astype(int)
events = pd.read_csv(B / "rmats/count_pass/fromGTF.MXE.txt", sep="\t").set_index("ID")
events = events.rename(columns={"1stExonStart_0base": "s1", "1stExonEnd": "e1", "2ndExonStart_0base": "s2", "2ndExonEnd": "e2"})
lr = lr.join(events[["geneSymbol", "chr", "strand", "s1", "e1", "s2", "e2", "upstreamEE", "downstreamES"]], on="id")
lr["lr_effect"] = lr.long_read_effect.astype(float)
lr = lr.loc[lr.strand.eq("+") & lr.lr_effect.ne(0)].sort_values("p_value").head(24)

junctions = pd.read_csv(B / "pseudobulk_junctions_junctions.tsv.gz", sep="\t", usecols=["chromosome", "start", "end"])
index = {(c, s, e): i for i, (c, s, e) in enumerate(zip(junctions.chromosome, junctions.start, junctions.end))}
counts = sparse.load_npz(B / "pseudobulk_junctions_counts.npz").tocsc()
samples = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t")
rows = []
for e in lr.itertuples(index=False):
    u, s1, e1, s2, e2, d = e.upstreamEE, e.s1, e.e1, e.s2, e.e2, e.downstreamES
    spec = {"U-1": (u + 1, s1), "1-D": (e1 + 1, d), "U-2": (u + 1, s2), "2-D": (e2 + 1, d), "1-2": (e1 + 1, s2)}
    out = dict(gene=e.geneSymbol, contrast=f"{e.level_a} vs {e.level_b}", p=e.p_value, short_effect=e.short_read_effect, lr_effect=e.lr_effect, agree=e.pooled_replicated)
    for side, level in (("a", e.level_a), ("b", e.level_b)):
        rows_i = np.flatnonzero(samples.cell_type.eq(level).to_numpy())
        n = {k: float(counts[rows_i, index[(e.chr, *v)]].sum()) if (e.chr, *v) in index else 0. for k, v in spec.items()}
        inc, exc = n["U-1"] + n["1-D"], n["U-2"] + n["2-D"]
        out.update({f"{k}_{side}": n[k] for k in spec})
        out[f"inc_share_{side}"] = inc / (inc + exc) if inc + exc else np.nan
        out[f"specific_inc_share_{side}"] = n["1-D"] / (n["1-D"] + n["U-2"]) if n["1-D"] + n["U-2"] else np.nan
        total = sum(n.values())
        out[f"both_share_{side}"] = n["1-2"] / total if total else np.nan
        out[f"specific_share_{side}"] = (n["1-D"] + n["U-2"]) / total if total else np.nan
    rows.append(out)
d = pd.DataFrame(rows); d.to_csv(sys.argv[1], sep="\t", index=False)
pd.set_option("display.width", 250)
print(d[["gene", "contrast", "short_effect", "lr_effect", "agree", "both_share_a", "specific_share_a", "inc_share_a", "inc_share_b", "specific_inc_share_a", "specific_inc_share_b"]].round(3).to_string(index=False))
print("median share of event junction reads on 1-2:", round(d[["both_share_a", "both_share_b"]].stack().median(), 3), "| on the form-specific junctions 1-D and U-2:", round(d[["specific_share_a", "specific_share_b"]].stack().median(), 3))
sign = lambda x: np.sign(x)
d["specific_dir"] = sign(d.specific_inc_share_b - d.specific_inc_share_a)
print("form-specific direction agrees with long reads:", int((d.specific_dir == sign(d.lr_effect)).sum()), "of", int(d.specific_dir.ne(0).sum() - d.specific_dir.isna().sum()), "| rMATS direction agrees:", int((sign(d.short_effect) == sign(d.lr_effect)).sum()), "of", len(d))
