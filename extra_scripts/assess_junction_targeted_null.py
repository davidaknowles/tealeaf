#!/usr/bin/env python3
"""Coverage-only null p-values of the small-p LeafCutter or scQuint tests of one contrast.

Features (LeafCutter clusters, scQuint intron groups) with p below --threshold in the contrast's real result are
redrawn --draws times: each sample of the contrast keeps its feature total, and its junction counts are multinomial
from the proportions pooled over the contrast's samples (no difference between the two cell types). Each draw's copy
of a feature is renamed (chromosome suffixed ~<draw>) so all draws are tested in one run of the method, with the
production settings (LeafCutter: leafcutter_ds.R with the real run's metadata and filters; scQuint:
run_junction_comparator.py scquint). Both methods test each feature on its own counts, so the subset run gives the
same test as a full run. Output: contrast_id, feature_id, p_value, p_null<k> (1 when a draw's copy is not tested).
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from tealeaf.sc.junction_benchmark import JunctionBundle  # noqa: E402

LEAFCUTTER_FILTERS = ["--min_samples_per_intron", "3", "--min_samples_per_group", "3", "--min_coverage"]


def redraw(matrix, rng):
    """matrix (junctions x samples) -> multinomial redraw per sample from the proportions pooled over the samples."""
    pooled = matrix.sum(1).astype(float)
    if pooled.sum() == 0:
        return matrix.copy()
    return rng.multinomial(matrix.sum(0).astype(np.int64), pooled / pooled.sum()).T


def leafcutter_null(args, contrast, selected, work, rng):
    names, rows = [], []
    with gzip.open(args.counts, "rt") as handle:
        header = handle.readline().rstrip("\n").split(" ")
        for line in handle:
            name, rest = line.split(" ", 1)
            chrom, _, _, cluster = name.split(":")
            if f"{chrom}:{cluster}" in selected:
                names.append(name)
                rows.append(np.array(rest.split(), dtype=np.int64))
    if not names:
        return {}
    counts = np.vstack(rows)
    columns = np.array([header.index(s) for s in pd.read_csv(args.metadata, sep="\t", header=None)[0]])
    key = np.array([":".join([n.split(":")[0], n.split(":")[3]]) for n in names])
    lines = [" ".join(header)]
    for draw in range(args.draws):
        null = counts.copy()
        for feature in np.unique(key):
            index = np.flatnonzero(key == feature)
            null[np.ix_(index, columns)] = redraw(counts[np.ix_(index, columns)], rng)
        for name, row in zip(names, null):
            chrom, start, end, cluster = name.split(":")
            lines.append(f"{chrom}~{draw}:{start}:{end}:{cluster} " + " ".join(map(str, row)))
    path = work / "null_perind_numers.counts.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n")
    script = (f"unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset \"$v\"; done; module purge; module load leafcutter/0.2.9; "
              f"Rscript {args.leafcutter_root}/scripts/leafcutter_ds.R --num_threads 1 {' '.join(LEAFCUTTER_FILTERS)} {args.leafcutter_min_coverage} --output_prefix {work}/null {path} {args.metadata}")
    subprocess.run(["bash", "-lc", script], check=False, capture_output=True)
    significance = work / "null_cluster_significance.txt"
    if not significance.exists():
        return {}
    table = pd.read_csv(significance, sep="\t")
    table = table.loc[pd.to_numeric(table.p, errors="coerce").notna()]
    out = {}
    for cluster, p in zip(table.cluster, table.p.astype(float)):
        chrom_draw, clu = cluster.split(":", 1)
        chrom, draw = chrom_draw.rsplit("~", 1)
        out[(f"{chrom}:{clu}", int(draw))] = p
    return out


def scquint_null(args, contrast, selected, work, rng):
    bundle = JunctionBundle.load(args.bundle)
    junctions = bundle.junctions
    use = junctions.intron_group.isin(selected) & (junctions.intron_group_size >= 2)
    index = np.flatnonzero(use.to_numpy())
    if not len(index):
        return {}
    samples = {s: i for i, s in enumerate(bundle.samples.sample_id)}
    columns = np.array([samples[s] for s in contrast["samples_a"] + contrast["samples_b"] if s in samples])
    counts = bundle.counts.tocsc()[:, index].toarray().T  # junctions x samples
    groups = junctions.intron_group.to_numpy()[index]
    parts, matrices = [], []
    for draw in range(args.draws):
        null = np.zeros_like(counts)
        for group in np.unique(groups):
            rows = np.flatnonzero(groups == group)
            null[np.ix_(rows, columns)] = redraw(counts[np.ix_(rows, columns)], rng)
        local = junctions.iloc[index].copy()
        local["chromosome"] = local.chromosome.astype(str) + f"~{draw}"
        local["intron_group"] = local.intron_group.astype(str) + f"~{draw}"
        parts.append(local)
        matrices.append(null)
    from scipy import sparse
    small = JunctionBundle(counts=sparse.csr_matrix(np.vstack(matrices).T), junctions=pd.concat(parts, ignore_index=True), samples=bundle.samples)
    small.save(work / "null")
    contrasts = work / "contrast.json"
    contrasts.write_text(json.dumps([contrast]))
    output = work / "null_scquint.tsv"
    subprocess.run([sys.executable, str(REPO / "extra_scripts/run_junction_comparator.py"), "scquint", "--bundle", str(work / "null"), "--contrasts", str(contrasts),
                    "--contrast-index", "0", "--scquint-root", str(args.scquint_root), "--jobs", "1", "--output", str(output)], check=True, env=os.environ | {"PYTHONPATH": str(REPO)})
    if output.stat().st_size == 0:
        return {}
    try:
        table = pd.read_csv(output, sep="\t")
    except pd.errors.EmptyDataError:
        return {}
    out = {}
    for feature, p in zip(table.feature_id, table.p_value):
        name, draw = str(feature).rsplit("~", 1)
        out[(name, int(draw))] = float(p) if np.isfinite(p) else 1.
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=("leafcutter", "scquint"), required=True)
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--contrast-index", type=int, required=True)
    parser.add_argument("--tests", type=Path, required=True, help="the contrast's normalized real result")
    parser.add_argument("--counts", type=Path, help="LeafCutter perind numerator counts")
    parser.add_argument("--metadata", type=Path, help="the contrast's LeafCutter metadata")
    parser.add_argument("--bundle", type=Path, help="junction bundle prefix (scQuint)")
    parser.add_argument("--leafcutter-root", type=Path, default=Path.home() / "projects/leafcutter")
    parser.add_argument("--scquint-root", type=Path, default=Path.home() / "projects/scquint")
    parser.add_argument("--leafcutter-min-coverage", type=int, default=20, help="leafcutter_ds --min_coverage, as in the real run")
    parser.add_argument("--threshold", type=float, default=.01)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    contrast = json.loads(args.contrasts.read_text())[args.contrast_index]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        real = pd.read_csv(args.tests, sep="\t")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        real = pd.DataFrame(columns=["feature_id", "p_value"])
    real = real.loc[real.p_value < args.threshold, ["feature_id", "p_value"]]
    columns = ["contrast_id", "feature_id", "p_value"] + [f"p_null{draw}" for draw in range(args.draws)]
    if real.empty:
        pd.DataFrame(columns=columns).to_csv(args.output, sep="\t", index=False)
        return
    rng = np.random.default_rng(args.contrast_index * 7919 + 1)
    keep = os.environ.get("TARGETED_NULL_WORK")
    with (tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", "/tmp")) if keep is None else open(os.devnull)) as tmp:
        work = Path(tmp) if keep is None else Path(keep)
        work.mkdir(parents=True, exist_ok=True)
        selected = set(real.feature_id.astype(str))
        nulls = (leafcutter_null if args.method == "leafcutter" else scquint_null)(args, contrast, selected, work, rng)
    for draw in range(args.draws):
        real[f"p_null{draw}"] = [nulls.get((str(feature), draw), 1.) for feature in real.feature_id]
    real.insert(0, "contrast_id", contrast["contrast_id"])
    real[columns].to_csv(args.output, sep="\t", index=False)
    print(f"{contrast['contrast_id']}: {len(real)} tests, median null p < .05: {int((real[[c for c in columns if c.startswith('p_null')]].median(axis=1) < .05).sum())}", flush=True)


if __name__ == "__main__":
    main()
