#!/usr/bin/env python3
"""Compare Bouchard CAVI, tilted VI, and CAVI-initialized tilted VI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from tealeaf.sc import ec_glmm, ec_glmm_full


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", type=int, default=6)
    parser.add_argument("--cavi-iterations", type=int, default=300)
    parser.add_argument("--initializer-iterations", type=int, default=25)
    parser.add_argument("--tilted-iterations", type=int, default=150)
    return parser.parse_args()


def simulated_data(seed, effect, *, observation_noise=False, isoforms=2):
    rng = np.random.default_rng(seed)
    mice = 10
    clusters = np.repeat(np.arange(mice), 2)
    condition = np.repeat(np.arange(mice) >= mice // 2, 2)
    cell_type = np.tile([0, 1], mice)
    design = np.column_stack((1 - cell_type, cell_type, condition))
    if isoforms == 2:
        mappings = (
            np.array([[1.0, 0.0], [0.35, 0.65], [0.0, 1.0]]),
            np.array([[1.0, 0.0], [0.70, 0.30], [0.0, 1.0]]),
        )
    else:
        mappings = (
            np.array([
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.55, 0.45, 0.0],
                [0.0, 0.35, 0.65],
            ]),
            np.array([
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.25, 0.75, 0.0],
                [0.0, 0.70, 0.30],
            ]),
        )
    dimension = isoforms - 1
    random_effect = rng.normal(0.0, 0.35, (mice, dimension))
    baseline = np.linspace(-0.3, 0.2, dimension)
    cell_shift = np.linspace(0.25, -0.15, dimension)
    logits = baseline[None, :] + cell_type[:, None] * cell_shift[None, :]
    logits[:, 0] += effect * condition
    logits += random_effect[clusters]
    if observation_noise:
        logits += rng.normal(0.0, 0.6, logits.shape)
    abundance = np.column_stack((np.exp(logits), np.ones(len(logits))))
    counts = []
    for mapping in mappings:
        probability = abundance @ mapping.T
        probability /= probability.sum(axis=1, keepdims=True)
        counts.append(
            np.asarray([rng.multinomial(180, row) for row in probability])
        )
    return ec_glmm.ECGLMMData(tuple(counts), mappings, design, clusters)


def timed_fit(callback):
    started = time.perf_counter()
    fit = callback()
    return fit, time.perf_counter() - started


def summarize(method, fit, seconds, tilted_objective=None):
    return {
        "method": method,
        "seconds": seconds,
        "converged": fit["converged"],
        "iterations": fit["iterations"],
        "objective": fit["objective"],
        "tilted_objective": (
            fit["objective"] if tilted_objective is None else tilted_objective
        ),
        "condition_coefficient": np.asarray(fit["coefficients"])[-1, 0],
        "random_effect_sd": fit["random_effect_sd"],
    }


def main():
    args = parse_args()
    rows = []
    scenarios = (
        ("random_intercept", False, 2),
        ("observation_noise", True, 3),
    )
    for scenario, observation_noise, isoforms in scenarios:
        for effect in (0.0, 0.9):
            for seed in range(args.seeds):
                data = simulated_data(
                    seed,
                    effect,
                    observation_noise=observation_noise,
                    isoforms=isoforms,
                )
                cavi, seconds = timed_fit(
                    lambda: ec_glmm_full.fit_bouchard_cavi(
                        data,
                        observation_noise=observation_noise,
                        max_iter=args.cavi_iterations,
                    )
                )
                cavi_tilted = ec_glmm_full.fit_variational(
                    data,
                    observation_noise=observation_noise,
                    initial=cavi["parameters"],
                    max_iter=0,
                )["objective"]
                row = summarize("cavi", cavi, seconds, cavi_tilted)
                row.update({"scenario": scenario, "seed": seed, "effect": effect})
                rows.append(row)

                tilted, seconds = timed_fit(
                    lambda: ec_glmm_full.fit_variational(
                        data,
                        objective="tilted",
                        observation_noise=observation_noise,
                        max_iter=args.tilted_iterations,
                    )
                )
                row = summarize("tilted", tilted, seconds)
                row.update({"scenario": scenario, "seed": seed, "effect": effect})
                rows.append(row)

                hybrid, seconds = timed_fit(
                    lambda: ec_glmm_full.fit_cavi_then_tilted(
                        data,
                        observation_noise=observation_noise,
                        cavi_max_iter=args.initializer_iterations,
                        max_iter=args.tilted_iterations,
                    )
                )
                row = summarize("cavi_then_tilted", hybrid, seconds)
                row.update({"scenario": scenario, "seed": seed, "effect": effect})
                rows.append(row)

    table = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    summary = (
        table.groupby(["scenario", "effect", "method"], sort=False)
        .agg(
            fits=("seed", "size"),
            convergence=("converged", "mean"),
            median_seconds=("seconds", "median"),
            median_iterations=("iterations", "median"),
            median_tilted_objective=("tilted_objective", "median"),
            coefficient_bias=(
                "condition_coefficient",
                lambda values: float(np.mean(values))
                - float(table.loc[values.index, "effect"].iloc[0]),
            ),
            coefficient_rmse=(
                "condition_coefficient",
                lambda values: float(
                    np.sqrt(
                        np.mean(
                            (
                                values
                                - table.loc[values.index, "effect"].to_numpy()
                            )
                            ** 2
                        )
                    )
                ),
            ),
        )
        .reset_index()
    )
    summary_path = args.output.with_suffix(".summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(json.dumps({
        "fits": len(table),
        "output": str(args.output),
        "summary": str(summary_path),
    }, indent=2))


if __name__ == "__main__":
    main()
