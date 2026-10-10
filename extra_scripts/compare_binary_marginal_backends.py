#!/usr/bin/env python3
"""Check the grid/AGHQ binary EC model against the reference quadrature backend."""

import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.path_marginal import binary_marginal_test, prepare_binary_ec_likelihood
from tealeaf.sc.path_marginal_grid import KAPPA_BOUNDS, binary_grid_test, prepare_grid_likelihood
from tealeaf.sc.path_marginal_quadrature import shared_prior_quadrature


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--test-ids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cached = pickle.load(args.candidate_cache.open("rb"))
    requested = list(pd.read_csv(args.test_ids, sep="\t").test_id)
    candidates = {candidate[0]: candidate for candidate in cached["candidates"]}
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    rows = []
    for test_id in requested:
        _, _, _, gene, transcripts, path_index, signatures, sample_rows, _, levels = candidates[test_id]
        local_metadata, _, labels = local_test_design(metadata, sample_rows, levels, "cell_type_pairwise")
        subjects = local_metadata.mouse.astype(str).to_numpy()
        base, _, _ = local_gene_data(tuple(matrix[sample_rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
        baseline = pooled_isoform_weights(base, max_iter=250)
        record = {"test_id": test_id}
        start = time.monotonic()
        grid = binary_grid_test(prepare_grid_likelihood(base, path_index, labels, subjects, baseline))
        record.update(grid_seconds=time.monotonic() - start, grid_statistic=grid["statistic"], grid_converged=grid["converged"], grid_kappa=grid["alternative_concentration"], grid_sigma=grid["alternative_subject_sd"], grid_error=grid["quadrature_error"])
        try:
            start = time.monotonic()
            reference = binary_marginal_test(shared_prior_quadrature(prepare_binary_ec_likelihood(base, path_index, labels, subjects, baseline=baseline)), subject_nodes=96, path_nodes=64)
            record.update(reference_seconds=time.monotonic() - start, reference_statistic=reference["statistic"], reference_converged=reference["converged"], reference_kappa=reference["alternative_concentration"], reference_sigma=reference["alternative_subject_sd"], reference_error=reference["quadrature_error"])
        except ValueError as exception:
            record["reference_error_message"] = str(exception)
        rows.append(record)
        print(record, flush=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False, na_rep="NA")
    print(f"kappa bounds for grid backend {KAPPA_BOUNDS}; reference allows 1e4", flush=True)


if __name__ == "__main__":
    main()
