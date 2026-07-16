#!/bin/bash
set -euo pipefail

SCRIPT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/run_microglia_single_cell_glm_gpu.sbatch
CACHE_SCRIPT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/run_microglia_fixed_ec_design.sbatch
METHODS=${METHODS:-"admm_factorized frank_wolfe factorized"}
DESIGNS=${DESIGNS:-"binary weighted"}
TARGETS=${TARGETS:-"theta"}
ALEVIN=${ALEVIN:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin}

CACHE_JOB=""
if [[ " ${DESIGNS} " == *" weighted "* ]]; then
  if [[ ! -f "${ALEVIN}/gene_eqclass_fixed_weights.npz" || \
        "${ALEVIN}/gene_eqclass_probs.tsv.gz" -nt "${ALEVIN}/gene_eqclass_fixed_weights.npz" ]]; then
    CACHE_JOB=$(sbatch --parsable "${CACHE_SCRIPT}")
  fi
fi

for design in ${DESIGNS}; do
  dependency=()
  if [[ "${design}" == "weighted" && -n "${CACHE_JOB}" ]]; then
    dependency=(--dependency="afterok:${CACHE_JOB}")
  fi
  for target in ${TARGETS}; do
    for method in ${METHODS}; do
      sbatch "${dependency[@]}" \
        --job-name="tealeaf_${design}_${method}_${target}" \
        --export=ALL,METHOD="${method}",EC_DESIGN="${design}",REGULARIZATION_TARGET="${target}" \
        "${SCRIPT}"
    done
  done
done
