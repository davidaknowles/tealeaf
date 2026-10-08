#!/usr/bin/env python3
"""Unconditional mixed-binomial null checks on the same frozen count draws."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy import stats

from extra_scripts.audit_conditional_read_null import null_read_panels
from tealeaf.sc.local_read_mixed import LocalReadMixed, local_read_mixed_test, MODEL_VERSION


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", type=int, required=True)
    parser.add_argument("--draws", type=int, default=64)
    parser.add_argument("--nodes", type=int, default=11)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or args.draws < 1 or not 0 <= args.scenario < 4:
        raise ValueError("new output and declared finite null panel required")
    rows = []
    for draw in range(args.draws):
        counts, conditioned = null_read_panels(args.scenario, draw)
        for law, observed in (("unconditional binomial", counts), ("conditional working model", conditioned)):
            start = time.monotonic()
            record = dict(scenario=args.scenario, draw=draw, law=law, p_value=1., F_reference_p_value=1., converged=False, error="")
            try:
                likelihood = LocalReadMixed(observed)
                result = local_read_mixed_test(likelihood, nodes=args.nodes)
                record.update(result)
                record["F_reference_p_value"] = float(stats.f.sf(result["statistic"], 1, result["n_subjects"] - 1)) if result["converged"] else 1.
            except (ValueError, np.linalg.LinAlgError) as exc:
                record["error"] = str(exc)
            record["runtime_seconds"] = time.monotonic() - start
            rows.append(record)
        if draw % 8 == 0:
            print(f"scenario {args.scenario}, draw {draw}, {len(rows)} requested trials complete", flush=True)
    table = pd.DataFrame(rows)
    summaries = []
    for law, local in table.groupby("law"):
        for column in ("p_value", "F_reference_p_value"):
            summaries.append(dict(law=law, p_value_column=column, requested=len(local), converged=int(local.converged.sum()), nominal_05=int(local[column].le(.05).sum()), nominal_01=int(local[column].le(.01).sum()), median_runtime_seconds=float(local.runtime_seconds.median())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(model_version=MODEL_VERSION, scenario=args.scenario, draws=args.draws, nodes=args.nodes, seed=[20261008, 92317, args.scenario], frozen_draw_source="same null_read_panels generator as the conditional candidate", model="primer fixed means, independent normal shared-subject baseline and subject effect, beta null zero; both variance regimes refitted including exact zero", primary_calibration_law="unconditional binomial draws match this latent model, baseline SD1 and slope SD0/.8", cross_model_law="conditional draws are a sensitivity, not this model's generating law; their fixed class margins can change intercept/slope dependence", failure_policy="all requested trials retained at p=1", scope="toy count-null and numerical/runtime pilot, not real EC validation or a complete split/LR endpoint result", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
