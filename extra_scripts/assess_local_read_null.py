#!/usr/bin/env python3
"""Collate whole finite-count toy null families without successful-fit filtering."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.conditional_read_odds import MODEL_VERSION as CONDITIONAL_VERSION
from tealeaf.sc.local_read_mixed import MODEL_VERSION as UNCONDITIONAL_VERSION


LAWS = ("unconditional binomial", "conditional working model")


def integration_recipe(manifest, model):
    """Normalize explicitly declared node fields from the two model drivers."""
    key = 'quadrature_nodes' if model == 'conditional' else 'nodes'
    nodes, adaptive = manifest[key], manifest.get('adaptive_integration', False)
    if not isinstance(nodes, int) or isinstance(nodes, bool) or nodes < 3 or not isinstance(adaptive, bool):
        raise ValueError('declared integer quadrature order and boolean integration policy required')
    return dict(nodes=nodes, adaptive_integration=adaptive)


def validate_trials(table, scenario, draws):
    """Failures retain p=1; require every original law/draw identity."""
    if table.duplicated(["law", "draw"]).any() or set(zip(table.law, table.draw)) != {(law, draw) for law in LAWS for draw in range(draws)} or not table.scenario.eq(scenario).all():
        raise ValueError("complete unique original law/draw family required")
    columns = ["p_value", "F_reference_p_value"]
    values = table[columns].to_numpy(dtype=float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any() or not table.converged.isin([True, False]).all():
        raise ValueError("finite declared p-values and explicit fitting availability required")
    if not table.loc[~table.converged, columns].eq(1.).all().all():
        raise ValueError("unavailable trials must remain p=1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conditional-root", type=Path, required=True)
    parser.add_argument("--unconditional-root", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=64)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or args.draws < 1:
        raise ValueError("new output and positive declared draw count required")
    summaries, failures, hashes, integration_recipes = [], [], {}, {}
    for model, root, version in (("conditional", args.conditional_root, CONDITIONAL_VERSION), ("unconditional", args.unconditional_root, UNCONDITIONAL_VERSION)):
        for scenario in range(4):
            folder = root / f"scenario{scenario}"
            paths = [folder / name for name in ("tests.tsv.gz", "manifest.json")]
            manifest = json.loads(paths[1].read_text())
            integration = integration_recipe(manifest, model)
            if not isinstance(integration['adaptive_integration'], bool) or model in integration_recipes and integration_recipes[model] != integration:
                raise ValueError('quadrature integration recipe must agree across all scenarios within a model')
            integration_recipes[model] = integration
            if manifest["model_version"] != version or manifest["scenario"] != scenario or manifest["draws"] != args.draws or manifest["seed"] != [20261008, 92317, scenario] or manifest["production_changes"] is not False:
                raise ValueError("matching frozen toy null recipe required")
            hashes.update({str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})
            table = pd.read_csv(paths[0], sep="\t")
            validate_trials(table, scenario, args.draws)
            for law, local in table.groupby("law"):
                primary = law == ("conditional working model" if model == "conditional" else "unconditional binomial")
                for column in ("p_value", "F_reference_p_value"):
                    summaries.append(dict(model=model, scenario=scenario, prevalence=.002 if scenario // 2 else .3, true_slope_sd=.8 if scenario % 2 else 0., law=law, primary_generating_law=primary, tail=column, requested=len(local), usable=int(local.converged.sum()), nominal_05=int(local[column].le(.05).sum()), nominal_01=int(local[column].le(.01).sum()), median_runtime_seconds=float(local.runtime_seconds.median())))
                unavailable = local.loc[~local.converged]
                for row in unavailable.itertuples(index=False):
                    failures.append(dict(model=model, scenario=scenario, law=law, draw=row.draw, exception=str(row.error) if pd.notna(row.error) else "", parameter_boundary=getattr(row, "parameter_boundary", np.nan), quadrature_error=getattr(row, "quadrature_error", np.nan), n_subjects=getattr(row, "n_subjects", np.nan)))
    args.output_dir.mkdir(parents=True)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(args.output_dir / "unavailable_trials.tsv", sep="\t", index=False)
    manifest = dict(input_hashes=hashes, integration_recipes=integration_recipes, requested_models=2, requested_scenarios=4, laws=LAWS, draws_per_law=args.draws, complete_trials=2 * 4 * 2 * args.draws, scope="toy local-marker null diagnostics only, no original EC nulls, real power, complete-gene FDR or either replication endpoint", caveats="The two laws differ in whether margins precede independent random effects. Matched unconditional draws share declared seeds/generator; no raw-count draw hashes were archived by the original drivers. Count-model assumptions, rare-path fitting failures and limited sample sizes prevent certification. F tails are sensitivities, not independently calibrated alternatives.", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    summary = pd.DataFrame(summaries)
    print(summary.loc[summary.primary_generating_law & summary['tail'].eq('p_value')].to_string(index=False), flush=True)
    print(pd.DataFrame(failures).groupby(["model", "scenario"]).size().to_string(), flush=True)


if __name__ == "__main__":
    main()
