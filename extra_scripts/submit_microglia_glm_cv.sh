#!/bin/bash
set -euo pipefail

ROOT=/gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts
tune_jobs=()
fit_jobs=()
for design in binary weighted; do
  for method in admm_factorized frank_wolfe_penalized; do
    tune=$(sbatch --parsable \
      --job-name="tealeaf_cv_${design}_${method}" \
      --export=ALL,DESIGN="${design}",METHOD="${method}" \
      "${ROOT}/run_microglia_glm_cv.sbatch")
    fit=$(sbatch --parsable --dependency="afterok:${tune}" \
      --job-name="tealeaf_cvfit_${design}_${method}" \
      --export=ALL,DESIGN="${design}",METHOD="${method}" \
      "${ROOT}/run_microglia_glm_cv_selected.sbatch")
    tune_jobs+=("${tune}")
    fit_jobs+=("${fit}")
    echo "Submitted tune ${tune} and full fit ${fit}: ${design} ${method}"
  done
done

dependencies=$(IFS=:; echo "${fit_jobs[*]}")
score=$(sbatch --parsable --dependency="afterok:${dependencies}" \
  "${ROOT}/run_microglia_score_single_cell_glm.sbatch")
echo "Submitted ${score}: dependent log-gene PCA scoring"
