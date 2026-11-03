"""Paired rMATS on real vs coverage-only null counts (all events, global pooled inclusion), finished contrasts."""
import glob, sys
from pathlib import Path
import pandas as pd
R = Path(sys.argv[1]); C = Path(sys.argv[2])
rows = []
real_all = {k: pd.read_csv(C / f"split_data_fold{k}_tests.tsv.gz", sep="\t", usecols=["contrast_id", "event_type", "p_value", "q_value"], low_memory=False) for k in (0,)}
for path in sorted(glob.glob(str(R / "split_null/results/*.tsv.gz"))):
    null = pd.read_csv(path, sep="\t", usecols=["contrast_id", "event_type", "p_value", "q_value"])
    cid = null.contrast_id.iloc[0]
    real = real_all[0].loc[real_all[0].contrast_id.eq(cid)]
    for label, t in (("real", real), ("null", null)):
        t = t[t.p_value.notna()]
        for et, d in t.groupby("event_type"):
            rows.append(dict(contrast=cid, input=label, event_type=et, tests=len(d), p05=(d.p_value < .05).mean(), bh05=int((d.q_value < .05).sum())))
r = pd.DataFrame(rows)
print(r.pivot_table(index=["contrast", "event_type"], columns="input", values=["bh05", "p05"]).round(4).to_string())
r.to_csv(sys.argv[3], sep="\t", index=False)
