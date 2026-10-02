"""Sensitivity test of subject-paired PSI rather than ILR in saved hybrid fits."""

import argparse
import json
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import signed_null_p_value
from tealeaf.sc.differential import paired_mean_test


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260927)
    args = parser.parse_args()
    records, nulls = [], []
    for shard in sorted(args.shards.glob("shard_*")):
        tests = pd.read_csv(shard / "paired_path.tsv", sep="\t")
        usage = pd.read_csv(shard / "path_usage.tsv.gz", sep="\t")
        grouped = dict(tuple(usage.groupby("test_id", sort=False)))
        for record in tests.to_dict("records"):
            local = grouped.get(record["test_id"])
            values = np.empty((0, 1))
            if local is not None:
                local = local.copy()
                local["cell_type"] = local.cell_type.astype(str).replace({"0": record["level_a"], "1": record["level_b"]})
                paired = local.pivot(index="subject", columns="cell_type", values="inclusion")
                if record["level_a"] in paired and record["level_b"] in paired:
                    values = (paired[record["level_b"]] - paired[record["level_a"]]).dropna().to_numpy()[:, None]
            record.update(paired_mean_test(values))
            mean = float(values.mean()) if len(values) else np.nan
            record.update(effect_size=mean, mean_difference_norm=abs(mean), test_coordinate="psi", path_pseudocount=record["report_pseudocount"])
            records.append(record)
            if record["converged"]:
                checksum = zlib.crc32(record["test_id"].encode())
                for replicate in range(args.null_replicates):
                    rng = np.random.default_rng(np.random.SeedSequence((args.seed, checksum, replicate)))
                    nulls.append({"test_id": record["test_id"], "block_id": record["block_id"], "replicate": replicate, "p_value": signed_null_p_value(values, None, rng, False, 0)})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame(nulls, columns=["test_id", "block_id", "replicate", "p_value"]).to_csv(args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(json.dumps([]) + "\n")


if __name__ == "__main__":
    main()
