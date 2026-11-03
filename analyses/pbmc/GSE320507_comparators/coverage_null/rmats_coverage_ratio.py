"""rMATS (PAIRADISE JCEC) call rate by coverage ratio between the levels (median per-sample IJC+SJC), per event type, fold 0."""
import json
import numpy as np, pandas as pd
B = "/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark"
contrasts = {c["contrast_id"]: c for c in json.load(open(f"{B}/rmats/split/fold0_contrasts.json"))}
t = pd.read_csv(f"{B}/comparator_suppa_rmats/rmats/split_data_fold0_tests.tsv.gz", sep="\t", usecols=["contrast_id", "event_type", "event_id", "p_value", "q_value", "effect_size"], low_memory=False)
t = t[t.p_value.notna()]
rows = []
for et in ("SE", "MXE", "A3SS", "A5SS"):
    raw = pd.read_csv(f"{B}/rmats/count_pass/JCEC.raw.input.{et}.txt", sep="\t")
    tot = np.array([list(map(int, a.split(","))) for a in raw.IJC_SAMPLE_1]) + np.array([list(map(int, s.split(","))) for s in raw.SJC_SAMPLE_1])
    pos = pd.Series(np.arange(len(raw)), index=raw.ID)
    d = t[t.event_type.eq(et)]
    ratio = np.empty(len(d))
    for cid, idx in d.groupby("contrast_id").groups.items():
        c = contrasts[cid]
        rr = pos[d.loc[idx, "event_id"]].to_numpy()
        a = np.median(tot[np.ix_(rr, c["indices_a"])], 1); b = np.median(tot[np.ix_(rr, c["indices_b"])], 1)
        ratio[d.index.get_indexer(idx)] = np.log2((a + 1) / (b + 1))
    d = d.assign(ratio=np.abs(ratio), call=d.q_value < .05)
    for lab, g in d.groupby(pd.cut(d.ratio, [-1, .5, 1, 2, 3, 100]), observed=True):
        rows.append(dict(event_type=et, abs_log2_coverage_ratio=str(lab), tests=len(g), call_rate=g.call.mean(), p05_rate=(g.p_value < .05).mean()))
r = pd.DataFrame(rows); print(r.to_string(index=False)); r.to_csv(".cache/pbmc/GSE320507/rmats_mxe_audit/call_rate_by_coverage_ratio.tsv", sep="\t", index=False)
g = pd.read_csv("analyses/pbmc/GSE320507_comparators/split/split_gene_tests.tsv.gz", sep="\t")
rm = set(g.loc[g.comparison.eq("rMATS (paired)") & g.method.eq("rMATS (paired)") & g.universe.eq("comparator") & (g.conjunction_q_value <= .05), "gene_id"])
f = pd.read_csv("analyses/pbmc/GSE320507_comparators/majiq_filtered_split/split_gene_tests.tsv.gz", sep="\t")
for v in ("all", "null_p_ge_05"):
    m = set(f.loc[f.comparison.eq(f"MAJIQ Heterogen [{v}]") & f.method.eq(f"MAJIQ Heterogen [{v}]") & f.universe.eq("comparator") & (f.conjunction_q_value <= .05), "gene_id"])
    print("MAJIQ", v, "replicated", len(m), "shared with rMATS", len(m & rm))
