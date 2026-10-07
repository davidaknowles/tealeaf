#!/usr/bin/env python3
"""Independent two-block challenge for fixed within-path isoform mixtures.

Four transcripts form the Cartesian product of binary A and B splice paths.
Only A is tested. Changing B alters within-A transcript shares, not A usage.
Both primers have known read-origin compatibility and opportunity weights.
Local tests retain A-origin reads only, either junctions or junctions+exon.
This is a controlled diagnostic, not a replacement real-data benchmark.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_marginal import binary_marginal_test, prepare_binary_ec_likelihood
from tealeaf.sc.path_marginal_quadrature import shared_prior_quadrature
from tealeaf.sc.path_simulation import simulate_independent_binary_blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=64)
    parser.add_argument("--subjects", type=int, default=12)
    parser.add_argument("--gene-depth", type=int, default=100)
    parser.add_argument("--a-effect", type=float, default=0.)
    parser.add_argument("--b-effect", type=float, default=0.)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=8)
    args = parser.parse_args()
    if args.draws < 1 or args.subjects < 4 or args.gene_depth < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("invalid size, depth or shard")
    output = []
    for draw in range(args.shard_index, args.draws, args.shard_count):
        rng = np.random.default_rng(7314159 + draw)
        generated, truth = simulate_independent_binary_blocks(rng, n_subjects=args.subjects, gene_depth=args.gene_depth, a_effect=args.a_effect, b_effect=args.b_effect)
        labels, subjects, paths, baseline = (truth[key] for key in ("labels", "subjects", "path_index", "baseline"))
        for strategy, selected in (("All gene reads, fixed within-path shares", np.arange(7)), ("Local junction and exon reads", np.arange(4)), ("Local junction reads", np.arange(3))):
            start = time.monotonic()
            try:
                data = ECGLMMData(tuple(values[:, selected] for values in generated.counts), tuple(mapping[selected] for mapping in generated.compatibility), generated.design, subjects)
                likelihood = shared_prior_quadrature(prepare_binary_ec_likelihood(data, paths, labels, subjects, baseline=baseline))
                fitted = binary_marginal_test(likelihood, subject_nodes=96, path_nodes=32)
                result = {key: fitted[key] for key in ("p_value", "statistic", "converged", "quadrature_error", "null_concentration", "alternative_concentration")}
                result["estimated_delta"] = fitted["standardized_means"][1, 0] - fitted["standardized_means"][0, 0]
                result["error"] = ""
            except (ValueError, np.linalg.LinAlgError) as exception:
                result = {"p_value": 1., "converged": False, "estimated_delta": np.nan, "error": str(exception)}
            output.append({"draw": draw, "strategy": strategy, "true_delta": truth["true_delta"], "runtime_seconds": time.monotonic() - start, **result})
            print(f"{draw}, {strategy}, p={result['p_value']:.4g}, converged={result['converged']}, delta={result['estimated_delta']:.4g}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "settings.json").write_text(json.dumps({**vars(args), "output_dir": str(args.output_dir), "true_subject_sd": .3, "biological_type_residuals": "none, same A usage within a subject when a-effect=0", "baseline": "label-blind mean of generating transcript weights, diagnostic oracle", "test_model": "same random-subject Beta EC LRT for all three read selections", "latent_model": "independent A/B path composition, four Cartesian-product transcripts", "production_changes": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
