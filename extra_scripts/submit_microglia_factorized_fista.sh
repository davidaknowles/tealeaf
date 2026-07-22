#!/bin/bash
set -euo pipefail

ROOT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts
fit_jobs=()
for design in binary weighted; do
  tune=$(sbatch --parsable --gres=gpu:b6k:1 \
    --job-name="tealeaf_cv_${design}_factor_fista" \
    --export=ALL,DESIGN="${design}",METHOD=factorized,CV_CELLS=0,MIN_CELL_UMIS=500,OUTPUT_TAG=_minumi500_all_fista,MAX_ITER=4096,MAX_RANK=512,EXACT_INNER_STEPS=32 \
    "${ROOT}/run_microglia_glm_cv.sbatch")
  fit=$(sbatch --parsable --gres=gpu:b6k:1 --dependency="afterok:${tune}" \
    --job-name="tealeaf_fit_${design}_factor_fista" \
    --export=ALL,DESIGN="${design}",METHOD=factorized,REPORT_TAG=_minumi500_all_fista,FIT_TAG=cv_fista,MIN_CELL_UMIS=500,FULL_MAX_ITER=8192,EXACT_INNER_STEPS=32 \
    "${ROOT}/run_microglia_glm_cv_selected.sbatch")
  fit_jobs+=("${fit}")
  echo "Submitted rank CV ${tune} and selected fit ${fit}: ${design}"
done

dependencies=$(IFS=:; echo "${fit_jobs[*]}")
score=$(sbatch --parsable --dependency="afterok:${dependencies}" \
  "${ROOT}/run_microglia_score_factorized_fista.sbatch")
echo "Submitted ${score}: dependent factorized FISTA scoring"
