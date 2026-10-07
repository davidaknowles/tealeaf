"""Validate and summarize the independent-block read-origin challenge."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


STRATEGIES = ("All gene reads, fixed within-path shares", "Local junction and exon reads", "Local junction reads")


def collate(cache):
    files = sorted(Path(cache).glob("shard_*/settings.json"))
    if not files:
        raise ValueError("no completed shards")
    settings = [json.loads(path.read_text()) for path in files]
    expected = settings[0]["shard_count"]
    if {path.parent.name for path in files} != {f"shard_{i}" for i in range(expected)}:
        raise ValueError("missing contiguous shards")
    common = [{key: value for key, value in item.items() if key not in ("output_dir", "shard_index")} for item in settings]
    if any(item != common[0] for item in common[1:]):
        raise ValueError("inconsistent shard settings")
    tables = []
    for path, setting in zip(files, settings):
        table = pd.read_csv(path.parent / "tests.tsv.gz", sep="\t")
        if not table.draw.mod(expected).eq(setting["shard_index"]).all():
            raise ValueError("incorrect shard assignment")
        tables.append(table)
    table = pd.concat(tables, ignore_index=True)
    strategies = common[0].get("strategies", STRATEGIES)
    keys = {(draw, strategy) for draw in range(common[0]["draws"]) for strategy in strategies}
    if table.duplicated(["draw", "strategy"]).any() or set(zip(table.draw, table.strategy)) != keys:
        raise ValueError("missing or duplicate draw/strategy keys")
    good = table.converged.astype(str).str.lower().eq("true")
    if not np.isfinite(table.loc[good, ["p_value", "estimated_delta"]]).all().all():
        raise ValueError("nonfinite successful result")
    table["converged"] = good
    table.loc[~good, "p_value"] = 1.
    table.loc[~good, "estimated_delta"] = np.nan
    table["effect_bias"] = table.estimated_delta - table.true_delta
    table["direction_agrees"] = np.where(good & table.true_delta.ne(0), np.sign(table.estimated_delta).eq(np.sign(table.true_delta)), np.nan)
    summary = table.groupby("strategy", sort=False).agg(n_requested=("p_value", "size"), n_converged=("converged", "sum"), native_reject_0_05=("p_value", lambda values: values.le(.05).sum()), native_reject_rate_0_05=("p_value", lambda values: values.le(.05).mean()), median_estimated_delta=("estimated_delta", "median"), mean_effect_bias=("effect_bias", "mean"), direction_agreement=("direction_agrees", "mean"), median_runtime_seconds=("runtime_seconds", "median")).reset_index()
    return table, summary, common[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    table, summary, settings = collate(args.cache)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps({"settings": settings, "scope": "Controlled independent-block misspecification challenge, not real-data split/LR performance", "failure_policy": "all requested trials retained, failed p=1 and effects unavailable", "production_changes": False}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
