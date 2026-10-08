#!/usr/bin/env python3
"""Separate conditional-model implementation checks from unconditional read nulls."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy import special, stats

from tealeaf.sc.conditional_read_odds import ConditionalReadOdds, conditional_read_odds_test, simulate_local_read_counts, simulate_conditional_read_counts, MODEL_VERSION
from tealeaf.sc.path_score_mixed import mixed_score_test


def gaussian_conditional_score(likelihood):
    """Count-table Gaussian sensitivity, not the existing global-EC method."""
    values = [likelihood.subject_evaluate(index, [0.]) for index in range(likelihood.n_informative_subjects)]
    score = np.array([row[1][0] for row in values])[:, None]
    information = np.array([row[2][0] for row in values])[:, None, None]
    return mixed_score_test(score, information, np.ones_like(information), reference_information=information, scalar_fast=True)


def null_read_panels(scenario, draw):
    """Bit-identical frozen toy draws for conditional/unconditional candidates."""
    if not 0 <= scenario < 4 or draw < 0:
        raise ValueError("declared scenario and nonnegative draw required")
    prevalence, sd = (.002 if scenario // 2 else .3), (.8 if scenario % 2 else 0.)
    rng = np.random.default_rng(np.random.SeedSequence((20261008, 92317, scenario, draw)))
    subjects = 20
    subject_capture = np.exp(rng.normal(0., .7, size=subjects))
    totals = np.maximum(1, np.rint(subject_capture[:, None, None] * np.array([[25, 100], [50, 20]])[None, :, :])).astype(int)
    baseline = special.logit(prevalence) + rng.normal(0., 1., size=(subjects, 1)) + np.array([1.5, -1.5])[None, :]
    effect = rng.normal(0., sd, size=subjects) if sd else np.zeros(subjects)
    counts = simulate_local_read_counts(totals, baseline, effect, rng)
    independent_effect = rng.normal(0., sd, size=subjects) if sd else np.zeros(subjects)
    conditioned = simulate_conditional_read_counts(counts, independent_effect, rng)
    return counts, conditioned


def run(args):
    if args.output_dir.exists() or args.draws < 1 or not 0 <= args.scenario < 4:
        raise ValueError("new output, positive draws and a declared scenario required")
    prevalence = .002 if args.scenario // 2 else .3
    sd = .8 if args.scenario % 2 else 0.
    rows = []
    for draw in range(args.draws):
        counts, conditioned = null_read_panels(args.scenario, draw)
        for law, observed in (("unconditional binomial", counts), ("conditional working model", conditioned)):
            start = time.monotonic()
            record = dict(scenario=args.scenario, prevalence=prevalence, true_subject_sd=sd, draw=draw, law=law, p_value=1., F_reference_p_value=1., score_p_value=1., statistic=0., converged=False, score_available=False, score_error="", error="")
            try:
                likelihood = ConditionalReadOdds(observed)
                record["n_informative_subjects"] = likelihood.n_informative_subjects
                if likelihood.n_informative_subjects >= 4:
                    try:
                        score = gaussian_conditional_score(likelihood)
                        record.update(score_p_value=score["p_value"], score_available=True)
                    except (ValueError, np.linalg.LinAlgError) as exc:
                        record["score_error"] = str(exc)
                result = conditional_read_odds_test(likelihood, nodes=args.nodes)
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
        for column in ("p_value", "F_reference_p_value", "score_p_value"):
            summaries.append(dict(law=law, p_value_column=column, requested=len(local), converged=int(local.converged.sum()), score_available=int(local.score_available.sum()), nominal_05=int(local[column].le(.05).sum()), nominal_01=int(local[column].le(.01).sum()), median_runtime_seconds=float(local.runtime_seconds.median())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(model_version=MODEL_VERSION, scenario=args.scenario, prevalence=prevalence, subject_sd=sd, draws=args.draws, subjects=20, primers=2, quadrature_nodes=args.nodes, seed=[20261008, 92317, args.scenario], failure_policy="requested trials remain p=1, never infer calibration from successful cases alone", target="mean local inclusion log-odds contrast zero, Gaussian subject effects shared across primers", laws="conditional draws assign a fresh independent effect to frozen margins; binomial draws condition only on type totals, not class margins", score="constant-log-odds-shape Gaussian conditional-count score sensitivity, NOT the existing global EC Tealeaf model", scope="finite unconditional/conditional local-marker toy count nulls, not real EC nulls, observed power, complete-gene FDR or either real replication endpoint", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", type=int, required=True)
    parser.add_argument("--draws", type=int, default=64)
    parser.add_argument("--nodes", type=int, default=21)
    parser.add_argument("--output-dir", type=Path, required=True)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
