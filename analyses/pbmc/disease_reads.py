"""Disease-group levels for GSE320507: relabel collection counts as <type>|<group> for unpaired within-type contrasts.

Groups are HC, CONV, LTCS (subject_folds.tsv condition), plus HC+CONV (pooled no long COVID) and a label-permutation
null (P12/P26: a random 12/26 split of the donors, seed 0) for calibration. Blocks, opportunities and filters are
shared with the source collection (symlinks). Writes <out>/{all,junction} and <out>/level_pairs.tsv.
Usage: disease_reads.py <dataset dir> <out dir> [<source collection name, default collection>]"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root, out = Path(sys.argv[1]), Path(sys.argv[2])
source = sys.argv[3] if len(sys.argv) > 3 else "collection"
folds = pd.read_csv(root / "subject_folds.tsv", sep="\t", dtype={"subject": str})
group = dict(zip(folds.subject, folds.condition))
pooled = {s: "HC+CONV" for s, g in group.items() if g in ("HC", "CONV")}
order = np.random.default_rng(0).permutation(sorted(group))
permuted = {s: "P12" if i < 12 else "P26" for i, s in enumerate(order)}
types = set()
for collection in ("all", "junction"):
    target = out / collection
    target.mkdir(parents=True, exist_ok=True)
    for name in ("blocks.json.gz", "opportunities.json.gz", "filters.json"):
        link = target / name
        if not link.exists():
            link.symlink_to((root / source / collection / name).resolve())
    counts = pd.read_csv(root / source / collection / "counts.tsv.gz", sep="\t", dtype={"subject": str})
    counts = counts.loc[counts.subject.isin(group)]
    types |= set(counts.cell_type)
    frames = []
    for labels in (group, pooled, permuted):
        local = counts.loc[counts.subject.isin(labels)].copy()
        local["cell_type"] = local.cell_type + "|" + local.subject.map(labels)
        frames.append(local)
    pd.concat(frames, ignore_index=True).to_csv(target / "counts.tsv.gz", sep="\t", index=False)
    print(collection, len(counts), "rows ->", sum(map(len, frames)))
contrasts = [("LTCS", "HC"), ("LTCS", "CONV"), ("CONV", "HC"), ("LTCS", "HC+CONV"), ("P26", "P12")]
pairs = pd.DataFrame([(f"{t}|{a}", f"{t}|{b}") for t in sorted(types) for a, b in contrasts], columns=["level_a", "level_b"])
pairs.to_csv(out / "level_pairs.tsv", sep="\t", index=False)
print(len(pairs), "level pairs")
