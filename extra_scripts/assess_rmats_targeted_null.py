#!/usr/bin/env python3
"""Coverage-only null p-values of the small-p paired rMATS (PAIRADISE, JCEC) tests of one contrast.

The contrast's per-sample inclusion/skipping counts come either from rMATS' prepare step on the all-sample count pass
(--count-pass) or from prepared per-contrast directories whose sample columns are concatenated in order
(--prepared, e.g. the two subject halves of a full-data contrast; PAIRADISE pairs sample i of the two groups and does
not depend on the order of the pairs). Events with p below --threshold in the real result are kept, and in each of
--draws draws every sample keeps its total I + S while I is redrawn Binomial(I + S, f0), f0 = sum I / sum (I + S)
pooled over the contrast's samples. PAIRADISE is run on the drawn counts with run_rmats_comparison.py's JCEC path.
Output: contrast_id, event_type, event_id, p_value, p_null<k> (1 when not tested).
"""

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_rmats_comparison import EVENT_TYPES, PREPARE, parse_result, run_command, run_paired_jcec_only  # noqa: E402


def read_counts(path):
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    split = lambda column: np.array([np.array(v.split(","), dtype=np.int64) if v else np.zeros(0, np.int64) for v in table[column]])
    return table, split("IJC_SAMPLE_1"), split("SJC_SAMPLE_1"), split("IJC_SAMPLE_2"), split("SJC_SAMPLE_2")


def join(values):
    return [",".join(map(str, row)) for row in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--contrast-index", type=int, required=True, help="position in --contrasts")
    parser.add_argument("--tests", type=Path, required=True, help="real paired rMATS tests (contrast_id, event_type, event_id, p_value)")
    parser.add_argument("--count-pass", type=Path)
    parser.add_argument("--prepared", type=Path, action="append", default=[])
    parser.add_argument("--threshold", type=float, default=.01)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--rscript", type=Path)
    parser.add_argument("--r-libs-user", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    contrast = json.loads(args.contrasts.read_text())[args.contrast_index]
    real = pd.concat(chunk.loc[chunk.contrast_id.eq(contrast["contrast_id"]) & (chunk.p_value < args.threshold)]
                     for chunk in pd.read_csv(args.tests, sep="\t", usecols=["contrast_id", "event_type", "event_id", "p_value"], chunksize=1_000_000))
    columns = ["contrast_id", "event_type", "event_id", "p_value"] + [f"p_null{draw}" for draw in range(args.draws)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if real.empty:
        pd.DataFrame(columns=columns).to_csv(args.output, sep="\t", index=False)
        return
    run_args = SimpleNamespace(rscript=args.rscript, r_libs_user=args.r_libs_user)
    rng = np.random.default_rng(args.contrast_index * 104729 + 17)
    nulls = {}
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", "/tmp")) as tmp:
        tmp = Path(tmp)
        sources = args.prepared
        if args.count_pass is not None:
            prepared = tmp / "prepared"
            run_command(["python", str(PREPARE), "--new-output-dir", prepared, "--old-output-dir", args.count_pass, "--group-1-indices", ",".join(map(str, contrast["indices_a"])), "--group-2-indices", ",".join(map(str, contrast["indices_b"]))], run_args)
            sources = [prepared]
        for event_type in EVENT_TYPES:
            ids = set(real.loc[real.event_type.eq(event_type), "event_id"].astype(int).astype(str))
            if not ids:
                continue
            parts = [read_counts(source / f"JCEC.raw.input.{event_type}.txt") for source in sources]
            table = parts[0][0]
            for other in parts[1:]:
                if not other[0].ID.equals(table.ID):
                    raise ValueError("prepared inputs have different events")
            inc_a, skip_a, inc_b, skip_b = (np.hstack([part[k] for part in parts]) for k in (1, 2, 3, 4))
            keep = table.ID.isin(ids).to_numpy()
            table, inc_a, skip_a, inc_b, skip_b = table.loc[keep].reset_index(drop=True), inc_a[keep], skip_a[keep], inc_b[keep], skip_b[keep]
            total_a, total_b = inc_a + skip_a, inc_b + skip_b
            f0 = (inc_a.sum(1) + inc_b.sum(1)) / np.maximum(total_a.sum(1) + total_b.sum(1), 1)
            gtf = pd.read_csv(sources[0] / f"fromGTF.{event_type}.txt", sep="\t", dtype=str)
            gtf = gtf.loc[gtf.ID.isin(ids)]
            for draw in range(args.draws):
                work = tmp / f"{event_type}_{draw}"
                work.mkdir()
                new_a = rng.binomial(total_a, f0[:, None])
                new_b = rng.binomial(total_b, f0[:, None])
                out = table.copy()
                out["IJC_SAMPLE_1"], out["SJC_SAMPLE_1"] = join(new_a), join(total_a - new_a)
                out["IJC_SAMPLE_2"], out["SJC_SAMPLE_2"] = join(new_b), join(total_b - new_b)
                out.to_csv(work / f"JCEC.raw.input.{event_type}.txt", sep="\t", index=False)
                gtf.to_csv(work / f"fromGTF.{event_type}.txt", sep="\t", index=False)
                run_paired_jcec_only(work, run_args, args.threads)
                result = parse_result(work, contrast)
                for event_id, p in zip(result.event_id.astype(str), result.p_value):
                    nulls[(event_type, event_id, draw)] = float(p) if np.isfinite(p) else 1.
                shutil.rmtree(work)
    real = real.assign(event_id=real.event_id.astype(int).astype(str))
    for draw in range(args.draws):
        real[f"p_null{draw}"] = [nulls.get((t, e, draw), 1.) for t, e in zip(real.event_type, real.event_id)]
    real[columns].to_csv(args.output, sep="\t", index=False)
    medians = real[[f"p_null{draw}" for draw in range(args.draws)]].median(axis=1)
    print(f"{contrast['contrast_id']}: {len(real)} tests, median null p < .05: {int((medians < .05).sum())}", flush=True)


if __name__ == "__main__":
    main()
