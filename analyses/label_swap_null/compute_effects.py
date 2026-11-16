"""Usage-shift effect sizes of the comparators' half-data tests in GSE320507, for the real labels and each swap draw.

LeafCutter (clusters) and scQuint (intron groups): the largest absolute change in junction share within the feature
between the two cell types, from junction counts pooled over each cell type's pseudobulks in the contrast (LeafCutter's
perind numerator counts; the scQuint junction bundle). rMATS: |IncLevelDifference| of its own output. Tealeaf effects
come from its fits (mean path usage difference) in perm_fdr_split.py.
Writes <output>: labeling, method, fold, pair_id, feature_id, effect.
Usage: compute_effects.py <output.parquet> <draw dir> [<draw dir> ...]"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
OUT, DRAWS = Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]]
key = lambda a, b: "||".join(sorted((a, b)))

# LeafCutter: perind numerators, junction rows chr:start:end:clu_x_strand, columns pseudobulk samples
lc = pd.read_csv(B / "leafcutter/clustering/benchmark_perind_numers.counts.gz", sep=" ", index_col=0, dtype={c: np.float32 for c in []})
lc_feature = [f"{n.split(':')[0]}:{n.split(':')[3]}" for n in lc.index]
lc_counts, lc_samples = lc.to_numpy(np.float32).T, list(lc.columns)  # samples x junctions
del lc
# scQuint: junction bundle, intron groups
sq_counts = sparse.load_npz(B / "pseudobulk_junctions_counts.npz").tocsr()
sq_samples = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t").sample_id.tolist()
sq_feature = pd.read_csv(B / "pseudobulk_junctions_junctions.tsv.gz", sep="\t", usecols=["intron_group"]).intron_group.astype(str).tolist()


def feature_shifts(counts, sample_index, features, samples_a, samples_b):
    """Largest |share_b - share_a| over the junctions of each feature, shares within the feature, counts pooled per side."""
    rows_a = [sample_index[s] for s in samples_a if s in sample_index]
    rows_b = [sample_index[s] for s in samples_b if s in sample_index]
    total = lambda rows: np.asarray(counts[rows].sum(axis=0)).ravel()
    frame = pd.DataFrame({"feature": features, "a": total(rows_a), "b": total(rows_b)})
    sums = frame.groupby("feature")[["a", "b"]].transform("sum")
    frame["shift"] = (frame.b / sums.b - frame.a / sums.a).abs()  # NaN when a side has no reads in the feature
    return frame.groupby("feature")["shift"].max()


rows = []
labelings = [("real", lambda k: B / f"reproducibility/fold{k}/contrasts.json")] + [(d.name, lambda k, d=d: d / f"fold{k}_junction_contrasts.json") for d in DRAWS]
lc_index, sq_index = {s: i for i, s in enumerate(lc_samples)}, {s: i for i, s in enumerate(sq_samples)}
for labeling, path in labelings:
    for k in (0, 1):
        for c in json.loads(path(k).read_text()):
            pair = key(c["level_a"], c["level_b"])
            for method, counts, index, features in (("LeafCutter", lc_counts, lc_index, lc_feature), ("scQuint", sq_counts, sq_index, sq_feature)):
                s = feature_shifts(counts, index, features, c["samples_a"], c["samples_b"])
                rows.append(pd.DataFrame({"labeling": labeling, "method": method, "fold": k, "pair_id": pair, "feature_id": s.index, "effect": s.to_numpy()}))
        print(labeling, k, "junction methods done", flush=True)
# rMATS
for labeling, source in [("real", None)] + [(d.name, d) for d in DRAWS]:
    for k in (0, 1):
        if source is None:
            frames = [pd.read_csv(B / f"comparator_suppa_rmats/rmats/split_data_fold{k}_tests.tsv.gz", sep="\t", usecols=["level_a", "level_b", "feature_id", "effect_size"])]
        else:
            frames = [pd.read_csv(p, sep="\t", usecols=["level_a", "level_b", "feature_id", "effect_size"]) for p in sorted((source / f"fold{k}/rmats").glob("*.tsv.gz"))]
        t = pd.concat(frames, ignore_index=True)
        rows.append(pd.DataFrame({"labeling": labeling, "method": "rMATS (paired)", "fold": k, "pair_id": [key(a, b) for a, b in zip(t.level_a, t.level_b)], "feature_id": t.feature_id, "effect": t.effect_size.abs()}))
        print(labeling, k, "rMATS done", flush=True)
out = pd.concat(rows, ignore_index=True)
out.to_parquet(OUT)
print(out.groupby(["labeling", "method"]).effect.describe().round(3).to_string())
