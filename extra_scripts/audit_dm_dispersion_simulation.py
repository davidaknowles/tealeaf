#!/usr/bin/env python3
"""Exact integer-count DM control, separating dispersion from EC approximation."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import softmax

from tealeaf.sc import differential
from tealeaf.sc.ec_block_glmm import blocked_multilevel_design
from tealeaf.sc.path_dispersion import corrected_dm_test


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=128)
    parser.add_argument("--subjects", type=int, default=12)
    parser.add_argument("--types", type=int, choices=(2, 5), default=2)
    parser.add_argument("--paths", type=int, choices=(2, 3), default=2)
    parser.add_argument("--concentration", type=float, default=20.)
    parser.add_argument("--effect", type=float, default=0.)
    parser.add_argument("--minimum-depth", type=int, default=200)
    parser.add_argument("--maximum-depth", type=int, default=2000)
    parser.add_argument("--mean-logit-offset", type=float, default=0.)
    parser.add_argument("--quantification-concentration", type=float, default=0., help="Optional uniform smoothing before fitting fractional effective counts, not extra reads.")
    parser.add_argument("--integerize", action="store_true", help="Diagnostic largest-remainder rounding after smoothing.")
    args = parser.parse_args()
    if args.draws < 1 or args.subjects < 4 or args.concentration <= 0 or not np.isfinite(args.concentration):
        parser.error("positive draws/concentration and at least four subjects required")
    if args.minimum_depth < 1 or args.maximum_depth < args.minimum_depth or not np.isfinite(args.mean_logit_offset) or not np.isfinite(args.quantification_concentration) or args.quantification_concentration < 0:
        parser.error("valid depth range, finite mean offset and nonnegative smoothing required")
    labels = np.tile(np.arange(args.types), args.subjects)
    subjects = np.repeat(np.arange(args.subjects), args.types)
    design, tested, _, _ = blocked_multilevel_design(labels, subjects)
    rng = np.random.default_rng(271828 + args.types * 100 + args.paths * 10000)
    subject_logits = rng.normal(0, .5, size=(args.subjects, args.paths - 1)) + args.mean_logit_offset
    type_effects = args.effect * np.linspace(-1, 1, args.types)[:, None] * np.linspace(1, -.5, args.paths - 1)[None, :]
    means = softmax(np.column_stack([subject_logits[subjects] + type_effects[labels], np.zeros(len(labels))]), axis=1)
    depths = np.tile(np.linspace(args.minimum_depth, args.maximum_depth, args.types).astype(int), args.subjects)
    rows = []
    for draw in range(args.draws):
        latent = np.asarray([rng.dirichlet(args.concentration * mean) for mean in means])
        counts = np.asarray([rng.multinomial(depth, probability) for depth, probability in zip(depths, latent)])
        if args.quantification_concentration:
            counts = depths[:, None] * (counts + args.quantification_concentration / args.paths) / (depths[:, None] + args.quantification_concentration)
        if args.integerize:
            counts = differential.integerize_compositional_counts(counts)
        results = {}
        for name, options in (("ML", None), ("Cox-Reid", {}), ("Known precision", {"concentration": args.concentration})):
            try:
                result = differential.dirichlet_multinomial_test(counts, design[:, :tested[0]], design, fix_null_concentration=False) if options is None else corrected_dm_test(counts, design[:, :tested[0]], design, **options)
                converged = bool(result["null_converged"] and result["alternative_converged"])
                results[name + " chi2"] = {**result, "p_value": result["p_value"] if converged else 1., "converged": converged, "error": ""}
                if options is not None:
                    results[name + " F"] = {**results[name + " chi2"], "p_value": result["f_p_value"] if converged else 1.}
            except (ValueError, np.linalg.LinAlgError) as exception:
                for tail in (("chi2",) if options is None else ("chi2", "F")):
                    results[name + " " + tail] = {"p_value": 1., "converged": False, "error": str(exception)}
        for strategy, result in results.items():
            rows.append({"draw": draw, "strategy": strategy, **{key: result.get(key, np.nan) for key in ("p_value", "converged", "error", "null_concentration", "alternative_concentration", "profile_boundary", "statistic", "degrees_of_freedom")}})
        if (draw + 1) % 16 == 0:
            print(f"{draw + 1}/{args.draws}", flush=True)
    output = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary = output.groupby("strategy").agg(n_tests=("p_value", "size"), n_converged=("converged", "sum"), reject_0_05=("p_value", lambda values: values.le(.05).mean()), reject_0_01=("p_value", lambda values: values.le(.01).mean()), precision_median=("alternative_concentration", "median")).reset_index()
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "settings.json").write_text(json.dumps({**vars(args), "output_dir": str(args.output_dir), "distribution": "exact subject-blocked integer DM counts, optionally transformed to uniformly smoothed fractional estimates before fitting", "same_model_as_fitted": args.quantification_concentration == 0, "null": args.effect == 0, "seed": 271828 + args.types * 100 + args.paths * 10000, "interpretation": "Native-tail and fractional-response diagnostic, not EC calibration or real-data endpoint"}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
