#!/usr/bin/env python3
"""Calibrate direct-count random-subject Beta inference on exact binary counts.

This uses the same low-depth/rare-path setting as the fixed-subject DM control.
Each trial redraws Gaussian subject effects and Beta subject/type compositions.
Primer counts share the composition. It is not a real-data endpoint assessment.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.special import expit

from tealeaf.sc.path_marginal import MODEL_VERSION, BinaryECPathLikelihood, binary_marginal_test


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=128)
    parser.add_argument("--subjects", type=int, default=12)
    parser.add_argument("--types", type=int, default=2)
    parser.add_argument("--concentration", type=float, default=20.)
    parser.add_argument("--subject-scale", type=float, default=.5)
    parser.add_argument("--mean-logit-offset", type=float, default=-2.2)
    parser.add_argument("--effect", type=float, default=0.)
    parser.add_argument("--minimum-depth", type=int, default=20)
    parser.add_argument("--maximum-depth", type=int, default=200)
    parser.add_argument("--subject-nodes", type=int, default=24)
    parser.add_argument("--path-nodes", type=int, default=16)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    args = parser.parse_args()
    if args.draws < 1 or args.subjects < 4 or args.types < 2 or args.minimum_depth < 1 or args.maximum_depth < args.minimum_depth or args.concentration <= 0 or args.subject_scale <= 0 or not 0 <= args.shard_index < args.shard_count:
        parser.error("invalid simulation size, depth, precision, subject scale or shard")
    labels = np.tile(np.arange(args.types), args.subjects)
    subjects = np.repeat(np.arange(args.subjects), args.types)
    depths = np.tile(np.linspace(args.minimum_depth, args.maximum_depth, args.types).astype(int), args.subjects)
    components = (np.array([[0., 1., 0.], [0., 0., 1.]]),) * 2
    output = []
    for draw in range(args.shard_index, args.draws, args.shard_count):
        # Seed per trial, independent of sharding and null-versus-effect choice.
        rng = np.random.default_rng(614159 + 100 * args.types + draw)
        mean = expit(args.mean_logit_offset + rng.normal(0, args.subject_scale, args.subjects)[subjects] + args.effect * np.linspace(-1, 1, args.types)[labels])
        latent = rng.beta(args.concentration * mean, args.concentration * (1 - mean))
        primer_depths = (depths // 2, depths - depths // 2)
        counts = tuple(np.column_stack([rng.binomial(depth, latent), np.zeros(len(labels))]) for depth in primer_depths)
        for values, depth in zip(counts, primer_depths):
            values[:, 1] = depth - values[:, 0]
        total = counts[0] + counts[1]
        proposals = depths[:, None] * (total + .125) / (depths[:, None] + .25)
        likelihood = BinaryECPathLikelihood(counts, components, subjects, labels, proposals, np.zeros(len(labels)))
        start = time.monotonic()
        try:
            fitted = binary_marginal_test(likelihood, subject_nodes=args.subject_nodes, path_nodes=args.path_nodes)
            result = {key: value for key, value in fitted.items() if key not in ("null_fit", "alternative_fit", "levels", "standardized_means")}
            result["effect"] = fitted["standardized_means"][-1, 0] - fitted["standardized_means"][0, 0]
            result["null_optimizer_success"] = bool(fitted["null_fit"].success)
            result["alternative_optimizer_success"] = bool(fitted["alternative_fit"].success)
            result["error"] = ""
        except (ValueError, np.linalg.LinAlgError) as exception:
            result = {"p_value": 1., "converged": False, "error": str(exception)}
        output.append({"draw": draw, "runtime_seconds": time.monotonic() - start, **result})
        print(f"{draw}, p={result['p_value']:.4g}, converged={result['converged']}, elapsed={output[-1]['runtime_seconds']:.1f}s", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    settings = {**vars(args), "output_dir": str(args.output_dir), "model_version": MODEL_VERSION, "quadrature_validation": "Both fitted objectives and their LR must agree on doubling both integration orders to absolute tolerance .001, no refitting at doubled order", "generating_model": "Gaussian random subject offset, Beta subject/type composition, independent primer binomials sharing composition", "unsupported": "multi-path blocks; no EC ambiguity in this synthetic control", "production_changes": False}
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")


if __name__ == "__main__":
    main()
