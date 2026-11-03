"""Calls of LeafCutter/scQuint on real vs coverage-only null inputs (full data, per-contrast BH as reported)."""
import sys
from pathlib import Path
import pandas as pd
P = Path(sys.argv[1]); rows = []


def read(f):
    try:
        return pd.read_csv(f, sep="\t", usecols=["contrast_id", "p_value", "q_value"])
    except pd.errors.EmptyDataError:
        print("empty", f)
        return None


for method, real, null in (("LeafCutter", P / "leafcutter/tests/normalized", P / "null_coverage/leafcutter/tests/normalized"),
                           ("scQuint", P / "scquint/tests", P / "null_coverage/scquint/tests")):
    for label, d in (("real", real), ("null", null)):
        files = sorted(d.glob("*.tsv"))
        t = pd.concat([x for x in map(read, files) if x is not None], ignore_index=True)
        t = t[t.p_value.notna()]
        rows.append(dict(method=method, input=label, contrasts=t.contrast_id.nunique(), tests=len(t), p05_rate=(t.p_value < .05).mean(),
                         p001_rate=(t.p_value < .001).mean(), bh05_calls=int((t.q_value < .05).sum())))
r = pd.DataFrame(rows); print(r.to_string(index=False)); r.to_csv(sys.argv[2], sep="\t", index=False)
