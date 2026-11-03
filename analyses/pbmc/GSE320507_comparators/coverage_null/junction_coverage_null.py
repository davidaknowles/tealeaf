"""Coverage-only null inputs for LeafCutter and scQuint on a junction benchmark.

Within each LeafCutter cluster (perind numerator counts) and each scQuint intron group (bundle intron_group, size >= 2),
every sample keeps its observed total and its junction counts are redrawn multinomially from the proportions pooled over
all samples. All cell types then share one usage, while expression (coverage) differences stay as observed.
Usage: junction_coverage_null.py <benchmark root> <output dir>
"""
import gzip, shutil, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tealeaf.sc.junction_benchmark import JunctionBundle  # noqa: E402

B, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261103)


def redraw(matrix):
    """matrix (junctions x samples) -> multinomial redraw per sample from pooled proportions."""
    totals = matrix.sum(0).astype(np.int64)
    pooled = matrix.sum(1).astype(float)
    if pooled.sum() == 0:
        return matrix
    draws = rng.multinomial(totals, pooled / pooled.sum())  # (samples, junctions)
    return draws.T


# LeafCutter perind numerators: header of sample names, rows "<chr:start:end:clu_N_strand> counts..."
src = B / "leafcutter/clustering/benchmark_perind_numers.counts.gz"
with gzip.open(src, "rt") as handle:
    header = handle.readline().rstrip("\n")
    names, rows = [], []
    for line in handle:
        fields = line.rstrip("\n").split(" ")
        names.append(fields[0])
        rows.append(np.array(fields[1:], dtype=np.int64))
counts = np.vstack(rows)
cluster = pd.Series(names).str.rsplit(":", n=1).str[1].to_numpy()
order = np.argsort(cluster, kind="stable")
bounds = np.flatnonzero(np.r_[True, cluster[order][1:] != cluster[order][:-1], True])
null = counts.copy()
for lo, hi in zip(bounds[:-1], bounds[1:]):
    rows_c = order[lo:hi]
    null[rows_c] = redraw(counts[rows_c])
with gzip.open(OUT / "leafcutter_perind_numers.counts.gz", "wt") as handle:
    handle.write(header + "\n")
    for name, row in zip(names, null):
        handle.write(name + " " + " ".join(map(str, row)) + "\n")
print("leafcutter clusters", len(bounds) - 1, "total", int(counts.sum()), "->", int(null.sum()), flush=True)

# scQuint: bundle counts (samples x junctions), redraw within intron groups of size >= 2
bundle = JunctionBundle.load(B / "pseudobulk_junctions")
X = bundle.counts.tocsc()
groups = bundle.junctions["intron_group"]
use = groups.notna() & (bundle.junctions["intron_group_size"] >= 2)
codes, uniq = pd.factorize(groups.where(use))
order = np.argsort(codes, kind="stable")
order = order[codes[order] >= 0]
c = codes[order]
bounds = np.flatnonzero(np.r_[True, c[1:] != c[:-1], True])
rows_i, cols_i, vals = [], [], []
for lo, hi in zip(bounds[:-1], bounds[1:]):
    cols = order[lo:hi]
    block = X[:, cols].toarray().T  # junctions x samples
    drawn = redraw(block)
    nz = np.nonzero(drawn)
    rows_i.append(nz[1]); cols_i.append(cols[nz[0]]); vals.append(drawn[nz])
keep = ~use.to_numpy()
rest = X[:, np.flatnonzero(keep)].tocoo()
rows_i.append(rest.row); cols_i.append(np.flatnonzero(keep)[rest.col]); vals.append(rest.data)
Y = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows_i), np.concatenate(cols_i))), shape=X.shape)
JunctionBundle(counts=Y, junctions=bundle.junctions, samples=bundle.samples).save(OUT / "pseudobulk_junctions")
print("scquint groups", len(bounds) - 1, "total", int(X.sum()), "->", int(Y.sum()), flush=True)
