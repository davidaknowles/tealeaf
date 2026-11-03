"""Audit of paired rMATS (PAIRADISE, JCEC) MXE calls in GSE320507 subject halves.

1. Calls by event type; MXE events and genes behind the calls.
2. Annotation of called MXE events vs all tested MXE events: transcripts containing exon 1, exon 2, both;
   whether the four junctions (upstream-exon1, exon1-downstream, upstream-exon2, exon2-downstream) are annotated introns.
3. Junction (JC) vs junction + exon body (JCEC) counts of called events, pooled over the contrast's samples per level:
   exon-body share of the inclusion/skipping counts, and whether the JC-only PSI difference has the sign of the
   reported effect.
Usage: rmats_mxe_audit.py <benchmark root> <gtf> <output dir>
"""
import gzip, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

B, GTF, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
R = B / "rmats"
C = B / "comparator_suppa_rmats/rmats"

# annotated exons and introns per transcript (GTF 1-based closed; rMATS starts are 0-based)
exons_tx = defaultdict(list)
with gzip.open(GTF, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "exon":
            continue
        tx = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        exons_tx[tx].append((f[0], int(f[3]) - 1, int(f[4])))
exon_tx = defaultdict(set)
introns = set()
for tx, ex in exons_tx.items():
    ex.sort(key=lambda e: e[1])
    for e in ex:
        exon_tx[e].add(tx)
    for a, b in zip(ex[:-1], ex[1:]):
        introns.add((a[0], a[2], b[1]))  # (chr, donor-side exon end, acceptor-side exon start 0-based)
print("transcripts", len(exons_tx), flush=True)

gtf_mxe = pd.read_csv(R / "count_pass/fromGTF.MXE.txt", sep="\t")
gtf_mxe.columns = [c.strip() for c in gtf_mxe.columns]
gtf_mxe["geneSymbol"] = gtf_mxe.geneSymbol.str.strip('"')


def annotate(row):
    chrom = row.chr
    e1 = (chrom, row["1stExonStart_0base"], row["1stExonEnd"])
    e2 = (chrom, row["2ndExonStart_0base"], row["2ndExonEnd"])
    t1, t2 = exon_tx.get(e1, set()), exon_tx.get(e2, set())
    j = [(chrom, row.upstreamEE, row["1stExonStart_0base"]), (chrom, row["1stExonEnd"], row.downstreamES),
         (chrom, row.upstreamEE, row["2ndExonStart_0base"]), (chrom, row["2ndExonEnd"], row.downstreamES)]
    return pd.Series({"tx_exon1": len(t1), "tx_exon2": len(t2), "tx_both": len(t1 & t2),
                      "annotated_junctions": sum(x in introns for x in j),
                      "exons_overlap": row["1stExonEnd"] > row["2ndExonStart_0base"],
                      "gap_exon1_exon2": row["2ndExonStart_0base"] - row["1stExonEnd"],
                      "len_exon1": row["1stExonEnd"] - row["1stExonStart_0base"], "len_exon2": row["2ndExonEnd"] - row["2ndExonStart_0base"]})


ann = pd.concat([gtf_mxe[["ID", "GeneID", "geneSymbol", "chr", "strand"]], gtf_mxe.apply(annotate, axis=1)], axis=1).rename(columns={"ID": "event_id"})


def counts(kind):
    t = pd.read_csv(R / f"count_pass/{kind}.raw.input.MXE.txt", sep="\t")
    inc = np.array([list(map(int, s.split(","))) for s in t.IJC_SAMPLE_1])
    skip = np.array([list(map(int, s.split(","))) for s in t.SJC_SAMPLE_1])
    return t.ID.to_numpy(), inc, skip, t.IncFormLen.to_numpy(float), t.SkipFormLen.to_numpy(float)


ids_jc, inc_jc, skip_jc, il_jc, sl_jc = counts("JC")
ids_ec, inc_ec, skip_ec, il_ec, sl_ec = counts("JCEC")
pos_jc = {e: i for i, e in enumerate(ids_jc)}
pos_ec = {e: i for i, e in enumerate(ids_ec)}

summary, gene_rows, call_rows = [], [], []
for k in (0, 1):
    contrasts = {c["contrast_id"]: c for c in json.loads((R / f"split/fold{k}_contrasts.json").read_text())}
    t = pd.read_csv(C / f"split_data_fold{k}_tests.tsv.gz", sep="\t", usecols=["contrast_id", "event_type", "event_id", "p_value", "q_value", "effect_size", "gene_name"], low_memory=False)
    t = t[t.p_value.notna()]
    for et, d in t.groupby("event_type"):
        summary.append(dict(fold=k, event_type=et, tests=len(d), calls=int((d.q_value < .05).sum()), call_rate=float((d.q_value < .05).mean()),
                            p05_rate=float((d.p_value < .05).mean())))
    calls = t[(t.q_value < .05) & t.event_type.eq("MXE")].copy()
    tested = t[t.event_type.eq("MXE")]
    calls["fold"] = k
    for row in calls.itertuples():
        c = contrasts[row.contrast_id]
        ia, ib = c["indices_a"], c["indices_b"]
        rec = dict(fold=k, contrast_id=row.contrast_id, event_id=row.event_id, effect=row.effect_size, q=row.q_value)
        for tag, pos, inc, skip, il, sl in (("jc", pos_jc, inc_jc, skip_jc, il_jc, sl_jc), ("jcec", pos_ec, inc_ec, skip_ec, il_ec, sl_ec)):
            i = pos[row.event_id]
            for lev, idx in (("a", ia), ("b", ib)):
                I, S = inc[i, idx].sum(), skip[i, idx].sum()
                rec[f"{tag}_inc_{lev}"], rec[f"{tag}_skip_{lev}"] = I, S
                rec[f"{tag}_psi_{lev}"] = (I / il[i]) / (I / il[i] + S / sl[i]) if I + S > 0 else np.nan
        call_rows.append(rec)
    # background annotation of tested vs called MXE events
    for label, d in (("tested", tested), ("called", calls)):
        a = ann.merge(d[["event_id"]].drop_duplicates(), on="event_id")
        summary.append(dict(fold=k, event_type=f"MXE {label} events (annotation)", tests=len(a),
                            calls=int(a.event_id.nunique()), call_rate=np.nan, p05_rate=np.nan,
                            frac_both_exons_one_transcript=float((a.tx_both > 0).mean()),
                            frac_exon1_unannotated=float((a.tx_exon1 == 0).mean()), frac_exon2_unannotated=float((a.tx_exon2 == 0).mean()),
                            frac_all4_junctions_annotated=float((a.annotated_junctions == 4).mean()),
                            frac_overlapping_exons=float(a.exons_overlap.mean()), median_gap=float(a.gap_exon1_exon2.median())))
    g = calls.merge(ann[["event_id", "geneSymbol"]], on="event_id").groupby("geneSymbol").agg(calls=("event_id", "size"), events=("event_id", "nunique"), contrasts=("contrast_id", "nunique")).sort_values("calls", ascending=False)
    g["fold"] = k
    gene_rows.append(g.reset_index())

summary = pd.DataFrame(summary)
genes = pd.concat(gene_rows, ignore_index=True)
calls = pd.DataFrame(call_rows).merge(ann, on="event_id", how="left")
calls["jc_effect"] = calls.jc_psi_a - calls.jc_psi_b
calls["jcec_effect"] = calls.jcec_psi_a - calls.jcec_psi_b
calls["body_share_inc"] = 1 - (calls.jc_inc_a + calls.jc_inc_b) / (calls.jcec_inc_a + calls.jcec_inc_b).replace(0, np.nan)
calls["body_share_skip"] = 1 - (calls.jc_skip_a + calls.jc_skip_b) / (calls.jcec_skip_a + calls.jcec_skip_b).replace(0, np.nan)
calls["jc_reads"] = calls[["jc_inc_a", "jc_inc_b", "jc_skip_a", "jc_skip_b"]].sum(1)
calls["jc_sign_agrees"] = np.sign(calls.jc_effect) == np.sign(calls.effect)
summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
genes.to_csv(OUT / "mxe_call_genes.tsv", sep="\t", index=False)
calls.to_csv(OUT / "mxe_calls.tsv.gz", sep="\t", index=False)
pd.set_option("display.width", 250)
print(summary.to_string(index=False))
print(genes.groupby("fold").head(15).to_string(index=False))
q = calls.groupby("fold").agg(calls=("event_id", "size"), events=("event_id", "nunique"),
                              median_body_share_inc=("body_share_inc", "median"), median_body_share_skip=("body_share_skip", "median"),
                              median_jc_reads=("jc_reads", "median"), frac_jc_reads_lt10=("jc_reads", lambda s: float((s < 10).mean())),
                              jc_sign_agrees=("jc_sign_agrees", "mean"), median_abs_effect=("effect", lambda s: float(s.abs().median())),
                              median_abs_jc_effect=("jc_effect", lambda s: float(s.abs().median())))
print(q.to_string())
