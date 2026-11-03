"""Coverage-only null count pass for paired rMATS (PAIRADISE).

For every event of every type, each sample keeps its total I_s + S_s (inclusion + skipping junction counts, JC and
JCEC separately) and I_s is redrawn Binomial(I_s + S_s, f0) with f0 = sum I / sum (I + S) pooled over all samples, so
all cell types share one inclusion level while coverage stays as observed. Other count-pass files are linked.
Usage: rmats_coverage_null.py <count pass dir> <output dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261103)
for path in src.iterdir():
    target = out / path.name
    if path.name.startswith(("JC.raw.input.", "JCEC.raw.input.")):
        t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        inc = np.array([np.array(v.split(","), dtype=np.int64) for v in t.IJC_SAMPLE_1])
        skip = np.array([np.array(v.split(","), dtype=np.int64) for v in t.SJC_SAMPLE_1])
        total = inc + skip
        f0 = inc.sum(1) / np.maximum(total.sum(1), 1)
        new_inc = rng.binomial(total, f0[:, None])
        t["IJC_SAMPLE_1"] = [",".join(map(str, row)) for row in new_inc]
        t["SJC_SAMPLE_1"] = [",".join(map(str, row)) for row in total - new_inc]
        t.to_csv(target, sep="\t", index=False)
        print(path.name, len(t), "events, totals kept", int(total.sum()), flush=True)
    elif not target.exists():
        target.symlink_to(path.resolve())
