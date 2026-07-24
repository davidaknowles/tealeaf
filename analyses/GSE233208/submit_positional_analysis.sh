#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"

ALEVIN_DIR="${DATA_ROOT}/processed/merged_alevin"
GLM_RUN="${RUN_ROOT}/glm"
mkdir -p "${GLM_RUN}/logs"
stamp=$(date +%Y%m%dT%H%M%S)
jobs_file="${GLM_RUN}/submitted_positional_jobs_${stamp}.tsv"
fits_file="${GLM_RUN}/submitted_positional_fits_${stamp}.tsv"
printf 'stage\tmethod\tjob_id\n' > "${jobs_file}"
: > "${fits_file}"

posbias=$(sbatch --parsable --array=0-39%10 \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  "${REPO_ROOT}/analyses/GSE233208/run_positional_bias.sbatch")
printf 'positional_design\tall\t%s\n' "${posbias}" >> "${jobs_file}"
design=$(sbatch --parsable --dependency="afterok:${posbias}" \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  "${REPO_ROOT}/analyses/GSE233208/build_positional_designs.sbatch")
printf 'positional_design_reduce\tall\t%s\n' "${design}" >> "${jobs_file}"
echo "positional quantification array=${posbias}; design reduction=${design}"

common="ALL,REPO_ROOT=${REPO_ROOT},ALEVIN_DIR=${ALEVIN_DIR},SALMON_REF=${SALMON_REF}/spliceu.fa,RUN_DIR=${GLM_RUN},CV_CELLS=0,MIN_CELL_UMIS=500,DESIGN=positional,ANALYSIS_TAG=paired_posbias,PRIMER_PAIRS=${ALEVIN_DIR}/primer_pairs.tsv,MIN_HALF_UMIS=500"
fit_jobs=()
for method in factorized admm_factorized frank_wolfe_penalized; do
  tune=$(sbatch --parsable --dependency="afterok:${design}" \
    --job-name="gse23_cv_pos_${method}" \
    -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
    --export="${common},METHOD=${method}" \
    "${REPO_ROOT}/extra_scripts/run_single_cell_glm_cv.sbatch")
  fit=$(sbatch --parsable --dependency="afterok:${tune}" \
    --job-name="gse23_fit_pos_${method}" \
    -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
    --export="${common},METHOD=${method}" \
    "${REPO_ROOT}/extra_scripts/run_single_cell_glm_selected.sbatch")
  fit_jobs+=("${fit}")
  printf 'cv\t%s\t%s\nfit\t%s\t%s\n' \
    "${method}" "${tune}" "${method}" "${fit}" >> "${jobs_file}"
  printf 'paired_posbias_positional_%s\t%s/paired_posbias_positional_%s_theta_\n' \
    "${method}" "${GLM_RUN}" "${method}" >> "${fits_file}"
  echo "${method}: CV=${tune} fit=${fit}"
done

score_jobs=()
for index in "${!fit_jobs[@]}"; do
  score=$(sbatch --parsable --dependency="afterok:${fit_jobs[index]}" \
    --job-name="gse23_score_pos_${index}" \
    -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
    --export="ALL,REPO_ROOT=${REPO_ROOT},FITS_FILE=${fits_file},FIT_INDEX=${index},LABELS=${DATA_ROOT}/metadata/reference_annotation.csv,GROUPS=${DATA_ROOT}/metadata/reference_donor.csv,TRANSCRIPT_TO_GENE=${SALMON_REF}/spliceu_t2g.tsv,SCORE_OUTPUT=${GLM_RUN}/label_scores" \
    "${REPO_ROOT}/extra_scripts/run_score_glm_fits.sbatch")
  score_jobs+=("${score}")
  printf 'score\t%s\t%s\n' "${index}" "${score}" >> "${jobs_file}"
done
score_dependency=$(IFS=:; echo "${score_jobs[*]}")
aggregate=$(sbatch --parsable --dependency="afterok:${score_dependency}" \
  --job-name=gse23_aggregate_pos --cpus-per-task=1 --mem=2G \
  --time=00:30:00 \
  -o "${GLM_RUN}/logs/%x-%j.out" -e "${GLM_RUN}/logs/%x-%j.err" \
  --wrap "bash -lc '${PYTHON_BIN} ${REPO_ROOT}/extra_scripts/aggregate_glm_scores.py --fits-file ${fits_file} --score-output ${GLM_RUN}/label_scores'")
printf 'aggregate\tall\t%s\n' "${aggregate}" >> "${jobs_file}"
echo "aggregate=${aggregate}"
echo "job manifest=${jobs_file}"
echo "fit manifest=${fits_file}"
