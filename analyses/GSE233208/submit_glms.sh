#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"

MERGE_JOB=${MERGE_JOB:?set MERGE_JOB to the successful/pending merge job id}
ALEVIN_DIR="${DATA_ROOT}/processed/merged_alevin"
GLM_RUN="${RUN_ROOT}/glm"
mkdir -p "${GLM_RUN}/logs"
stamp=$(date +%Y%m%dT%H%M%S)
jobs_file="${GLM_RUN}/submitted_jobs_${stamp}.tsv"
fits_file="${GLM_RUN}/submitted_fits_${stamp}.tsv"
printf 'stage\ttag\tdesign\tmethod\tjob_id\n' > "${jobs_file}"
: > "${fits_file}"

common="ALL,REPO_ROOT=${REPO_ROOT},ALEVIN_DIR=${ALEVIN_DIR},SALMON_REF=${SALMON_REF}/spliceu.fa,RUN_DIR=${GLM_RUN},CV_CELLS=0,MIN_CELL_UMIS=500"
validation=$(sbatch --parsable --dependency="afterok:${MERGE_JOB}" \
  --job-name=gse233208_validate \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  "${REPO_ROOT}/analyses/GSE233208/validate_merged.sbatch")
printf 'validate\tall\tall\tall\t%s\n' "${validation}" >> "${jobs_file}"
echo "merged data validation=${validation}"
cache=$(sbatch --parsable --dependency="afterok:${validation}" \
  --job-name=gse23_cache_designs --cpus-per-task=8 --mem=192G \
  --time=2-00:00:00 \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  --wrap "bash -lc 'export PYTHONPATH=${REPO_ROOT} && ${PYTHON_BIN} ${REPO_ROOT}/extra_scripts/cache_alevin_glm_designs.py --alevin-dir ${ALEVIN_DIR} --salmon-ref ${SALMON_REF}/spliceu.fa --primer-pairs ${ALEVIN_DIR}/primer_pairs.tsv --min-half-umis 500'")
echo "fixed weighted design cache=${cache}"
printf 'cache\tall\tall\tall\t%s\n' "${cache}" >> "${jobs_file}"
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
    printf 'cv\t%s\t%s\t%s\t%s\n' \
      "${tag}" "${design}" "${method}" "${tune}" >> "${jobs_file}"
    printf 'fit\t%s\t%s\t%s\t%s\n' \
      "${tag}" "${design}" "${method}" "${fit}" >> "${jobs_file}"
    printf '%s_%s_%s\t%s/%s_%s_%s_theta_\n' \
      "${tag}" "${design}" "${method}" "${GLM_RUN}" \
      "${tag}" "${design}" "${method}" >> "${fits_file}"
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
  printf 'cv\t%s\t%s\tfactorized\t%s\n' \
    "${tag}" "${design}" "${tune}" >> "${jobs_file}"
  printf 'fit\t%s\t%s\tfactorized\t%s\n' \
    "${tag}" "${design}" "${fit}" >> "${jobs_file}"
  printf '%s_%s_factorized\t%s/%s_%s_factorized_theta_\n' \
    "${tag}" "${design}" "${GLM_RUN}" "${tag}" "${design}" >> "${fits_file}"
  echo "${tag} ${design} factorized: CV=${tune} fit=${fit}"
done

score_jobs=()
for index in "${!fit_jobs[@]}"; do
  score=$(sbatch --parsable --dependency="afterok:${fit_jobs[index]}" \
    --job-name="gse23_score_${index}" \
    -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
    --export="ALL,REPO_ROOT=${REPO_ROOT},FITS_FILE=${fits_file},FIT_INDEX=${index},LABELS=${DATA_ROOT}/metadata/reference_annotation.csv,GROUPS=${DATA_ROOT}/metadata/reference_donor.csv,TRANSCRIPT_TO_GENE=${SALMON_REF}/spliceu_t2g.tsv,SCORE_OUTPUT=${GLM_RUN}/label_scores" \
    "${REPO_ROOT}/extra_scripts/run_score_glm_fits.sbatch")
  score_jobs+=("${score}")
  printf 'score\tall\tall\tfit_%s\t%s\n' \
    "${index}" "${score}" >> "${jobs_file}"
  echo "score fit ${index}=${score}"
done
score_dependency=$(IFS=:; echo "${score_jobs[*]}")
aggregate=$(sbatch --parsable --dependency="afterok:${score_dependency}" \
  --job-name=gse23_aggregate_scores --cpus-per-task=1 --mem=2G \
  --time=00:30:00 \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  --wrap "bash -lc '${PYTHON_BIN} ${REPO_ROOT}/extra_scripts/aggregate_glm_scores.py --fits-file ${fits_file} --score-output ${GLM_RUN}/label_scores'")
printf 'aggregate_scores\tall\tall\tall\t%s\n' \
  "${aggregate}" >> "${jobs_file}"
echo "aggregate scores=${aggregate}"
echo "job manifest=${jobs_file}"
echo "fit manifest=${fits_file}"
