"""Compare scalar/generic REML on the SAME archived count-null scores."""

import argparse
import json
from pathlib import Path
import zlib

import numpy as np
import pandas as pd

from extra_scripts.benchmark_event_score_kernels import compare_scalar


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    parts = [pd.read_csv(args.input_root / f"shard_{index}/subject_null_diagnostics.tsv.gz", sep="\t") for index in range(16)]
    table = pd.concat(parts, ignore_index=True)
    if table.duplicated(["test_id", "draw", "subject"]).any():
        raise ValueError("duplicate subject scores")
    rows = []
    for (test_id, draw), records in table.groupby(["test_id", "draw"], sort=False):
        scores = records.score.to_numpy()[:, None]
        info = records.information.to_numpy()[:, None, None]
        shapes = records.biological_shape.to_numpy()[:, None, None]
        reference = records.reference_information.to_numpy()[:, None, None]
        rng = np.random.default_rng(381924 + zlib.crc32(test_id.encode()) + 1721 * draw)
        rows.extend({"test_id": test_id, "count_draw": draw, **row} for row in compare_scalar(scores, info, shapes, reference, test_id, rng, 32))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.output_dir / "same_score_equivalence.tsv.gz", sep="\t", index=False)
    complete = frame.p_value_matches.notna()
    summary = dict(trials=len(frame), available=int(complete.sum()), same_availability=bool(frame.same_availability.all()), interpretation="same score, information, reference and biological shape arrays; no count or null refitting between implementations")
    for key in ("p_value", "statistic", "mean_difference", "mean_covariance", "restricted_objective", "biological_variance"):
        summary[key + "_all_match"] = bool(frame.loc[complete, key + "_matches"].all())
        summary[key + "_max_absolute_difference"] = float(frame[key + "_max_absolute_difference"].max())
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
