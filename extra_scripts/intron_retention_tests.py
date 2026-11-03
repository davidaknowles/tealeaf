#!/usr/bin/env python3
"""Differential intron retention between cell types: random-subject Beta-binomial test with a gene offset.

Counts come from intron_retention_counts.py. Per clean intron, row (subject, cell type), retained r is the
intron-body molecules (sense strand; poly(dT) molecules within --mask bases upstream of a sense A-run are
dropped, random-hexamer molecules are kept) and spliced n the molecules with the intron's junction. The offset
is the logit of (R + .5) / (N + .5), with R and N the same row's retained and spliced molecules over the gene's
other clean introns (intron_retention.gene_baseline_offset). Rows need r + n >= 3 and R + N >= 3, and a test
needs at least 4 subjects observed in both cell types. Every pair of --cell-types is tested
(intron_retention.retention_test).

--cohort is a fold of --folds or "full" (every subject of --folds). With --null-fraction f, that random fraction of
tests is fitted with each subject's cell-type labels swapped with probability 1/2, for the null calibration.
Writes tests_<cohort>[_null][_shard<k>].tsv.gz (concatenate shards before evaluation) with the likelihood-ratio statistic, the label effect (logit), the number of
subjects and the label-blind pooled retained and spliced molecules.
"""
import argparse
import itertools
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.intron_retention import gene_baseline_offset, retention_test

BINS = np.array([0, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 5000, np.inf])
MIN_SAMPLE, MIN_SUBJECTS = 3, 4


def load(counts_dir, subjects, cell_types, mask, retained="body"):
    keep_bin = np.searchsorted(BINS, mask, side="left")
    introns = pd.read_csv(counts_dir / "introns.tsv.gz", sep="\t")
    introns = introns[introns.clean].drop_duplicates(["chromosome", "start", "end"])
    parts = []
    for path in sorted(counts_dir.glob("counts_*.tsv.gz")):
        counts = pd.read_csv(path, sep="\t")
        counts = counts[counts.sense & counts.intron.isin(introns.intron) & counts.cell_type.isin(cell_types) & counts.subject.isin(subjects)]
        counts = counts[counts.kind.ne("body") | counts.primer.eq("random hexamer") | counts.bin.ge(keep_bin)]
        parts.append(counts.groupby(["intron", "subject", "cell_type", "kind"])["count"].sum())
    t = pd.concat(parts).groupby(level=[0, 1, 2, 3]).sum().unstack("kind", fill_value=0).reset_index()
    for kind in ("body", "left", "right", "spliced"):
        if kind not in t:
            t[kind] = 0
    t = t.merge(introns[["intron", "gene_id"]], on="intron")
    t["r"] = t.body + (t.left + t.right if retained == "body+boundary" else 0)
    t["n"] = t.spliced
    gene = t.groupby(["gene_id", "subject", "cell_type"])[["r", "n"]].transform("sum")
    t["R"], t["N"] = gene.r - t.r, gene.n - t.n
    t = t[(t.r + t.n).ge(MIN_SAMPLE) & (t.R + t.N).ge(MIN_SAMPLE)]
    t["offset"] = gene_baseline_offset(t.R, t.N)
    rows = []
    for a, b in itertools.combinations(sorted(cell_types), 2):
        m = t[t.cell_type.isin((a, b))]
        both = m.groupby(["intron", "subject"]).cell_type.transform("nunique").eq(2)
        m = m[both].assign(pair=f"{a}||{b}", label=lambda f: f.cell_type.eq(b).astype(int))
        m = m[m.groupby("intron").subject.transform("nunique").ge(MIN_SUBJECTS)]
        rows.append(m[["intron", "gene_id", "pair", "subject", "label", "r", "n", "offset"]])
    return pd.concat(rows, ignore_index=True)


def fit(group):
    (intron, pair), frame = group
    record = {"intron": intron, "pair": pair, "gene_id": frame.gene_id.iat[0], "n_subjects": frame.subject.nunique(), "retained": frame.r.sum(), "spliced": frame.n.sum(), "statistic": 0., "effect": np.nan, "p_chi2": 1., "converged": False, "kappa": np.nan, "sigma": np.nan}
    try:
        res = retention_test(frame.r.to_numpy(), frame.n.to_numpy(), frame.subject.to_numpy(), frame.label.to_numpy(), offset=frame.offset.to_numpy())
        record.update(statistic=res["statistic"], effect=float(res["coefficients"][1]), p_chi2=res["p_value"], converged=res["converged"], kappa=res["alternative_concentration"], sigma=res["alternative_subject_sd"])
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exception:
        record["error"] = str(exception)
    return record


def swapped(groups, fraction, seed):
    rng = np.random.default_rng(seed)
    result = []
    for i in np.flatnonzero(rng.random(len(groups)) < fraction):
        key, frame = groups[i]
        flip = dict(zip(frame.subject.unique(), rng.random(frame.subject.nunique()) < .5))
        result.append((key, frame.assign(label=frame.label ^ frame.subject.map(flip).astype(int))))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--folds", type=Path, required=True, help="tsv with subject and fold")
    parser.add_argument("--cohort", required=True, help="fold value or full")
    parser.add_argument("--cell-types", required=True, help="comma-separated cell types")
    parser.add_argument("--mask", type=float, default=1000)
    parser.add_argument("--retained", choices=("body", "body+boundary"), default="body")
    parser.add_argument("--null-fraction", type=float, default=0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shard", default="0/1", help="k/n fits every n-th (intron, pair) test from k; output gets a _shard<k> suffix when n > 1")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    folds = pd.read_csv(args.folds, sep="\t").set_index("subject").fold.astype(str)
    subjects = set(folds.index) if args.cohort == "full" else set(folds.index[folds.eq(args.cohort)])
    rows = load(args.counts_dir, subjects, args.cell_types.split(","), args.mask, args.retained)
    groups = list(rows.groupby(["intron", "pair"], sort=False))
    if args.null_fraction:
        groups = swapped(groups, args.null_fraction, seed=sum(map(ord, args.cohort)))
    k, n = map(int, args.shard.split("/"))
    groups = groups[k::n]
    print(args.cohort, "tests", len(groups), flush=True)
    with Pool(args.workers) as pool:
        records = pool.map(fit, groups, chunksize=200)
    pd.DataFrame(records).to_csv(args.output_dir / f"tests_{args.cohort}{'_null' if args.null_fraction else ''}{f'_shard{k}' if n > 1 else ''}.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
