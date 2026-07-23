#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"

MERGE_JOB=${MERGE_JOB:?set MERGE_JOB to the successful/pending merge job id}
ALEVIN_DIR="${DATA_ROOT}/processed/merged_alevin"
GLM_RUN="${RUN_ROOT}/glm"
mkdir -p "${GLM_RUN}/logs"

common="ALL,REPO_ROOT=${REPO_ROOT},ALEVIN_DIR=${ALEVIN_DIR},SALMON_REF=${SALMON_REF}/spliceu.fa,RUN_DIR=${GLM_RUN},CV_CELLS=0,MIN_CELL_UMIS=500"
cache=$(sbatch --parsable --dependency="afterok:${MERGE_JOB}" \
  --job-name=gse23_cache_designs --cpus-per-task=8 --mem=192G \
  --time=2-00:00:00 \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  --wrap "bash -lc 'module load python/3.10.8-GCCcore-12.2.0 && export PYTHONPATH=${REPO_ROOT} && python ${REPO_ROOT}/extra_scripts/cache_alevin_glm_designs.py --alevin-dir ${ALEVIN_DIR} --salmon-ref ${SALMON_REF}/spliceu.fa --primer-pairs ${ALEVIN_DIR}/primer_pairs.tsv --min-half-umis 500'")
echo "fixed weighted design cache=${cache}"
fit_jobs=()
for design in binary weighted; do
  for method in factorized admm_factorized frank_wolfe_penalized; do
    tag=standard
    tune=$(sbatch --parsable --dependency="afterok:${cache}" \
      --job-name="gse23_cv_${design}_${method}" \
      -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
      --export="${common},DESIGN=${design},METHOD=${method},ANALYSIS_TAG=${tag}" \
      "${REPO_ROOT}/extra_scripts/run_single_cell_glm_cv.sbatch")
    fit=$(sbatch --parsable --dependency="afterok:${tune}" \
      --job-name="gse23_fit_${design}_${method}" \
      -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
      --export="${common},DESIGN=${design},METHOD=${method},ANALYSIS_TAG=${tag}" \
      "${REPO_ROOT}/extra_scripts/run_single_cell_glm_selected.sbatch")
    fit_jobs+=("${fit}")
    echo "${tag} ${design} ${method}: CV=${tune} fit=${fit}"
  done

  tag=paired
  primer="${ALEVIN_DIR}/primer_pairs.tsv"
  tune=$(sbatch --parsable --dependency="afterok:${cache}" \
    --job-name="gse23_cv_${design}_paired" \
    -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
    --export="${common},DESIGN=${design},METHOD=factorized,ANALYSIS_TAG=${tag},PRIMER_PAIRS=${primer},MIN_HALF_UMIS=500" \
    "${REPO_ROOT}/extra_scripts/run_single_cell_glm_cv.sbatch")
  fit=$(sbatch --parsable --dependency="afterok:${tune}" \
    --job-name="gse23_fit_${design}_paired" \
    -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
    --export="${common},DESIGN=${design},METHOD=factorized,ANALYSIS_TAG=${tag},PRIMER_PAIRS=${primer},MIN_HALF_UMIS=500" \
    "${REPO_ROOT}/extra_scripts/run_single_cell_glm_selected.sbatch")
  fit_jobs+=("${fit}")
  echo "${tag} ${design} factorized: CV=${tune} fit=${fit}"
done
