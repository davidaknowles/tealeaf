#!/bin/bash
set -euo pipefail

SCRIPT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/run_microglia_single_cell_glm_gpu.sbatch
METHODS=${METHODS:-"admm_factorized frank_wolfe factorized"}

for method in ${METHODS}; do
  sbatch --export=ALL,METHOD="${method}" "${SCRIPT}"
done
