#!/usr/bin/env python3
"""Compare block-local class counts of two collections over the same blocks.

Used to compare the native target-locus aligner (native_collect.py, classified
into a STAR-derived block set) with the STARsolo-based collection
(collect_local_path_reads.py). Reports total molecules, the overlap
sum(min)/sum of the counts per (block, subject, cell type, primer, class),
the correlation of per-block totals, and the total-variation distance of each
block's class distribution (pooled over subjects and cell types).
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def summary(first, second):
    keys = ["path_key", "subject", "cell_type", "primer", "mask"]
    a = pd.read_csv(first / "counts.tsv.gz", sep="\t").groupby(keys)["count"].sum()
    b = pd.read_csv(second / "counts.tsv.gz", sep="\t").groupby(keys)["count"].sum()
    joined = pd.concat([a.rename("first"), b.rename("second")], axis=1).fillna(0)
    blocks = joined.groupby(level="path_key").sum()
    blocks = blocks[(blocks["first"] > 0) | (blocks["second"] > 0)]
    log_ratio = np.log2((blocks["second"] + 1) / (blocks["first"] + 1))
    classes = joined.groupby(level=["path_key", "mask"]).sum()
    shares = classes / classes.groupby(level="path_key").transform("sum")
    distance = (shares["first"] - shares["second"]).abs().groupby(level="path_key").sum() / 2
    weight = classes.groupby(level="path_key")["first"].sum()
    return {
        "molecules_first": float(joined["first"].sum()), "molecules_second": float(joined["second"].sum()),
        "overlap_of_first": float(np.minimum(joined["first"], joined["second"]).sum() / joined["first"].sum()),
        "overlap_of_second": float(np.minimum(joined["first"], joined["second"]).sum() / joined["second"].sum()),
        "blocks": int(len(blocks)), "block_total_pearson_log": float(np.corrcoef(np.log1p(blocks["first"]), np.log1p(blocks["second"]))[0, 1]),
        "block_log2_ratio_quantiles": dict(zip(("5%", "25%", "50%", "75%", "95%"), np.quantile(log_ratio, [.05, .25, .5, .75, .95]).round(3).tolist())),
        "class_distribution_tv_quantiles": dict(zip(("50%", "75%", "90%", "99%"), np.quantile(distance.dropna(), [.5, .75, .9, .99]).round(4).tolist())),
        "class_distribution_tv_weighted_mean": float((distance * weight).sum() / weight[distance.index].sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", nargs=3, action="append", metavar=("NAME", "FIRST", "SECOND"), required=True, help="name, reference collection (STAR), compared collection (native)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {name: summary(Path(first), Path(second)) for name, first, second in args.pair}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
