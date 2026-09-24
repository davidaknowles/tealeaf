#!/usr/bin/env python3
"""Run the rMATS event model for every paired benchmark contrast.

The all-sample rMATS count pass is intentionally separate from the per-contrast
statistical passes.  This avoids rescanning the BAMs while retaining rMATS'
native JCEC counts and unpaired beta-binomial test.  ``--output`` contains one
long table suitable for the shared split-reproducibility benchmark.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd


EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")
RMATS = Path("/nfs/sw/rmats/rmats-4.1.2/run_rmats")
PREPARE = Path("/nfs/sw/rmats/rmats-4.1.2/rMATS_P/prepare_stat_inputs.py")
OLD_OUTPUT = Path("/tmp/rmats_all/output")
LIB_DIR = Path("/tmp/rmats_lib")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contrasts", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path, default=Path("/tmp/rmats_stats"))
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--threads-per-job", type=int, default=2)
    parser.add_argument("--old-output", type=Path, default=OLD_OUTPUT)
    return parser.parse_args()


def command_environment():
    environment = os.environ.copy()
    old = environment.get("LD_LIBRARY_PATH", "")
    environment["LD_LIBRARY_PATH"] = f"{LIB_DIR}:{old}" if old else str(LIB_DIR)
    return environment


def run_command(command, cwd=None):
    # run_rmats has a fixed Python 3.9 environment supplied by the module.
    shell_command = "module load rmats/4.1.2 >/dev/null 2>&1; " + " ".join(map(str, command))
    subprocess.run(["bash", "-lc", shell_command], cwd=cwd, env=command_environment(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_result(path, contrast):
    rows = []
    for event_type in EVENT_TYPES:
        result_path = path / f"{event_type}.MATS.JCEC.txt"
        if not result_path.exists():
            continue
        table = pd.read_csv(result_path, sep="\t", low_memory=False, na_values=["NA", "nan", "NaN"])
        # The count-file ID is a second, duplicate ID column and is irrelevant.
        table = table.loc[:, ~table.columns.duplicated()]
        required = {"ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference"}
        if not required.issubset(table.columns):
            continue
        table = table.loc[:, ["ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference"]].copy()
        table["GeneID"] = table["GeneID"].astype(str).str.strip('"').str.split(";").str[0]
        table["geneSymbol"] = table["geneSymbol"].astype(str).str.strip('"').str.split(";").str[0]
        table["event_type"] = event_type
        table["event_id"] = table["ID"].astype(str)
        table["contrast_id"] = contrast["contrast_id"]
        table["effect"] = "cell_type"
        table["stratum"] = contrast.get("stratum", "all")
        table["level_a"] = contrast["level_a"]
        table["level_b"] = contrast["level_b"]
        table["method"] = "rMATS"
        table["feature_id"] = "rMATS:" + event_type + ":" + table["event_id"]
        table["p_value"] = pd.to_numeric(table["PValue"], errors="coerce")
        table["q_value"] = pd.to_numeric(table["FDR"], errors="coerce")
        table["effect_size"] = pd.to_numeric(table["IncLevelDifference"], errors="coerce")
        table["significant"] = table["q_value"].lt(0.05)
        table["criterion"] = "rMATS JCEC FDR < 0.05"
        table["gene_id"] = table["GeneID"].replace("nan", np.nan)
        table["gene_name"] = table["geneSymbol"].replace("nan", np.nan)
        rows.append(table[["method", "contrast_id", "effect", "stratum", "level_a", "level_b", "feature_id", "event_type", "event_id", "p_value", "q_value", "effect_size", "significant", "criterion", "gene_id", "gene_name"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def process_one(item):
    fold, contrast, old_output, work_root, threads = item
    index = int(contrast["index"])
    work = work_root / f"fold{fold}_{index}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    try:
        run_command(["python", str(PREPARE), "--new-output-dir", work, "--old-output-dir", old_output, "--group-1-indices", ",".join(map(str, contrast["indices_a"])), "--group-2-indices", ",".join(map(str, contrast["indices_b"]))])
        run_command([RMATS, "--od", work, "--tmp", work / "tmp", "--task", "stat", "--nthread", str(threads)], cwd=work)
        result = parse_result(work, contrast)
        result["fold"] = fold
        result.to_csv(work / "parsed.tsv", sep="\t", index=False)
        return fold, index, len(result), str(work / "parsed.tsv"), None
    except Exception as error:  # pragma: no cover - reported by parent
        return fold, index, 0, None, repr(error)


def main():
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    items = []
    for fold, path in enumerate(args.contrasts):
        contrasts = json.loads(path.read_text())
        for contrast in contrasts:
            items.append((fold, contrast, args.old_output, args.work_dir, args.threads_per_job))
    parsed = []
    failures = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(process_one, item) for item in items]
        for number, future in enumerate(as_completed(futures), 1):
            fold, index, count, path, error = future.result()
            print(f"{number}/{len(items)} fold={fold} contrast={index} rows={count}", flush=True)
            if error:
                failures.append({"fold": fold, "index": index, "error": error})
            elif count:
                parsed.append(path)
    if failures:
        (args.output.parent / "rmats_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
        raise RuntimeError(f"{len(failures)} rMATS contrasts failed; see rmats_failures.json")
    combined = pd.concat((pd.read_csv(path, sep="\t") for path in sorted(parsed)), ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(combined):,} rows to {args.output}")


if __name__ == "__main__":
    main()
