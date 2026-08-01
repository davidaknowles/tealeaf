#!/usr/bin/env python3
"""Simulate ambiguous EC counts and benchmark EC-GLMM inference methods."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from extra_scripts import run_ec_glmm
from tealeaf.sc import ec_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--effect-sizes", nargs="+", type=float, default=(0.0, 0.5))
    parser.add_argument("--mice", type=int, default=12)
    parser.add_argument("--total", type=int, default=150)
    parser.add_argument("--random-effect-sd", type=float, default=0.4)
    parser.add_argument("--concentration", type=float, default=30.0)
    parser.add_argument("--max-iter", type=int, default=150)
    parser.add_argument("--vi-samples", type=int, default=128)
    parser.add_argument("--renyi-alpha", type=float, default=0.5)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def simulate(args, family, effect, seed):
    rng = np.random.default_rng(seed)
    clusters = np.repeat(np.arange(args.mice), 2)
    condition = np.repeat(np.arange(args.mice) >= args.mice // 2, 2)
    cell_type = np.tile([0, 1], args.mice)
    null_design = np.c_[1 - cell_type, cell_type]
    alternative_design = np.c_[null_design, condition]
    random_effect = rng.normal(0.0, args.random_effect_sd, args.mice)
    logits = -0.35 * null_design[:, 0] + 0.25 * null_design[:, 1]
    logits += float(effect) * condition + random_effect[clusters]
    abundance = np.c_[np.exp(logits), np.ones(len(logits))]
    mappings = (
        np.array([[1.0, 0.0], [0.35, 0.65], [0.0, 1.0]]),
        np.array([[1.0, 0.0], [0.70, 0.30], [0.0, 1.0]]),
    )
    counts = []
    for mapping in mappings:
        probabilities = abundance @ mapping.T
        probabilities /= probabilities.sum(axis=1, keepdims=True)
        if family == "dirichlet_multinomial":
            probabilities = np.asarray([
                rng.dirichlet(args.concentration * row) for row in probabilities
            ])
        counts.append(np.asarray([
            rng.multinomial(args.total, row) for row in probabilities
        ]))
    null = ec_glmm.ECGLMMData(tuple(counts), mappings, null_design, clusters)
    alternative = ec_glmm.ECGLMMData(
        tuple(counts), mappings, alternative_design, clusters
    )
    return null, alternative


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    method_groups = {
        "multinomial": (
            "laplace_multinomial",
            "tilted_elbo",
            "renyi_multinomial",
        ),
        "dirichlet_multinomial": (
            "laplace_dirichlet_multinomial",
            "elbo_dirichlet_multinomial",
            "renyi_dirichlet_multinomial",
        ),
    }
    rows = []
    failures = []
    replicate_ids = range(args.shard_index, args.replicates, args.shard_count)
    for replicate in replicate_ids:
        for effect_number, effect in enumerate(args.effect_sizes):
            for family_number, (family, methods) in enumerate(method_groups.items()):
                seed = (
                    args.seed
                    + 100_003 * replicate
                    + 1_009 * effect_number
                    + 53 * family_number
                )
                null_data, alternative_data = simulate(args, family, effect, seed)
                starts = {}
                for method in methods:
                    started = time.perf_counter()
                    try:
                        null, alternative = run_ec_glmm.fit_nested(
                            method, null_data, alternative_data, args, starts
                        )
                        starts[method] = (null, alternative)
                        gain = (
                            2.0 * (null["objective"] - alternative["objective"])
                            if method.startswith("laplace")
                            else 2.0 * (alternative["objective"] - null["objective"])
                        )
                        rows.append({
                            "replicate": replicate,
                            "family": family,
                            "method": method,
                            "effect": effect,
                            "estimate": alternative["coefficients"][-1, 0],
                            "evidence_gain": gain,
                            "null_converged": null["converged"],
                            "alternative_converged": alternative["converged"],
                            "null_random_effect_sd": null["random_effect_sd"],
                            "alternative_random_effect_sd": alternative[
                                "random_effect_sd"
                            ],
                            "null_concentration": null["concentration"],
                            "alternative_concentration": alternative[
                                "concentration"
                            ],
                            "null_importance_ess": null.get("importance_ess", np.nan),
                            "alternative_importance_ess": alternative.get(
                                "importance_ess", np.nan
                            ),
                            "seconds": time.perf_counter() - started,
                        })
                    except Exception as exc:
                        failures.append({
                            "replicate": replicate,
                            "family": family,
                            "method": method,
                            "effect": effect,
                            "error": repr(exc),
                        })
        print(json.dumps({
            "replicate": replicate,
            "rows": len(rows),
            "failures": len(failures),
        }), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    args.output.with_suffix(".failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
