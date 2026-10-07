"""Diagnose within-path nuisance sharing in the existing path estimators.

All strategies use the same unmoderated paired mean test of fitted ILRs,
not the production family calibration. No residual biological variation
is added. These controlled tests diagnose estimators, not family FDR.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from tealeaf.sc import differential
from tealeaf.sc.path_simulation import simulate_independent_binary_blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=64)
    parser.add_argument("--a-effect", type=float, default=0.)
    parser.add_argument("--b-effect", type=float, default=0.)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=8)
    args = parser.parse_args()
    if args.draws < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("invalid number of trials or shard")
    strategies = (("Fixed within-path shares, A1", differential.fit_path_perturbation, 1.), ("Fixed within-path shares, A32", differential.fit_path_perturbation, 32.), ("Free within-path shares, A1", differential.fit_free_isoform_paths, 1.))
    rows = []
    for draw in range(args.shard_index, args.draws, args.shard_count):
        data, truth = simulate_independent_binary_blocks(np.random.default_rng(7314159 + draw), a_effect=args.a_effect, b_effect=args.b_effect)
        for strategy, fitter, concentration in strategies:
            start = time.monotonic()
            try:
                fits = [fitter(tuple(count[row] for count in data.counts), data.compatibility, truth["baseline"], truth["path_index"], path_pseudocount=concentration, path_pseudocount_scaling="total") for row in range(len(truth["labels"]))]
                if not all(fit.converged and fit.covariance.identifiable for fit in fits):
                    raise ValueError("an observation fit failed or is unidentifiable")
                ilrs = np.array([fit.path_logratios for fit in fits]).reshape(12, 2, -1)
                usage = np.array([fit.path_proportions for fit in fits]).reshape(12, 2, -1)
                test = differential.paired_mean_test(ilrs[:, 1] - ilrs[:, 0])
                if not test["converged"]:
                    raise ValueError("paired mean test failed")
                result = {key: test[key] for key in ("p_value", "statistic")}
                result.update(converged=True, estimated_delta=float(np.mean(usage[:, 1, 0] - usage[:, 0, 0])), error="")
            except (ValueError, np.linalg.LinAlgError) as exception:
                result = {"p_value": 1., "converged": False, "estimated_delta": np.nan, "error": str(exception)}
            rows.append({"draw": draw, "strategy": strategy, "true_delta": truth["true_delta"], "runtime_seconds": time.monotonic() - start, **result})
            print(f"{draw}, {strategy}, p={result['p_value']:.4g}, converged={result['converged']}, delta={result['estimated_delta']:.4g}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    settings = {**vars(args), "output_dir": str(args.output_dir), "strategies": [item[0] for item in strategies], "subjects": 12, "gene_depth": 100, "subject_sd": .3, "test": "unmoderated paired mean test of ILRs, no production calibration", "baseline": "label-blind generating pooled baseline, diagnostic oracle", "biological_type_residuals": "none", "production_changes": False}
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")


if __name__ == "__main__":
    main()
