#!/bin/bash
set -euo pipefail

SCRIPT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/run_microglia_single_cell_glm_gpu.sbatch
CACHE_SCRIPT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/run_microglia_fixed_ec_design.sbatch
METHODS=${METHODS:-"admm_factorized frank_wolfe factorized"}
DESIGNS=${DESIGNS:-"binary weighted"}
TARGETS=${TARGETS:-"theta"}

CACHE_JOB=""
if [[ " ${DESIGNS} " == *" weighted "* ]]; then
  CACHE_JOB=$(sbatch --parsable "${CACHE_SCRIPT}")
fi

for design in ${DESIGNS}; do
  dependency=()
  if [[ "${design}" == "weighted" ]]; then
    dependency=(--dependency="afterok:${CACHE_JOB}")
  fi
  for target in ${TARGETS}; do
    for method in ${METHODS}; do
      sbatch "${dependency[@]}" \
        --export=ALL,METHOD="${method}",EC_DESIGN="${design}",REGULARIZATION_TARGET="${target}" \
        "${SCRIPT}"
    done
  done
done
