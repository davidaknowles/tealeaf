#!/bin/bash
set -euo pipefail

source "$(dirname "$0")/config.env"

PARTITION=${PARTITION:-cpu}
TIME_LIMIT=${TIME_LIMIT:-02:00:00}
SHARD_LIMIT=${SHARD_LIMIT:-48}
NULL_REPLICATES=${NULL_REPLICATES:-32}
STRENGTHS=${STRENGTHS:-2;4;8;16;32;64;128;256}
LOG_ROOT="${RUN_ROOT}/junction_benchmark/reproducibility/logs"
RUN_FULL_OMNIBUS=${RUN_FULL_OMNIBUS:-1}
RUN_FULL_CONDITION=${RUN_FULL_CONDITION:-1}
RUN_FULL_PAIRWISE=${RUN_FULL_PAIRWISE:-1}

submit_family() {
  local family=$1
  local candidates=$2
  local moderation_flag=$3
  local result_root="${RUN_ROOT}/differential/local_path_full_testing_eb_${family}"
  local merge_dependencies=()
  local strengths
  local strength
  IFS=';' read -r -a strengths <<< "${STRENGTHS}"
  mkdir -p "${result_root}"
  for strength in "${strengths[@]}"; do
    local output_root="${result_root}/a${strength}_shards"
    local merged="${result_root}/a${strength}"
    local fit_job
    local merge_job
    fit_job=$(sbatch \
      --parsable \
      -p "${PARTITION}" \
      --time="${TIME_LIMIT}" \
      --array="0-95%${SHARD_LIMIT}" \
      --output="${LOG_ROOT}/full_testing_eb_${family}_a${strength}.%A_%a.out" \
      --error="${LOG_ROOT}/full_testing_eb_${family}_a${strength}.%A_%a.err" \
      --export="ALL,CONFIG=${REPO_ROOT}/analyses/microglia_less/config.env,CANDIDATE_CACHE=${candidates},OUTPUT_ROOT=${output_root},PATH_PSEUDOCOUNT=${strength},PATH_PRIOR_CENTER=uniform,PATH_PSEUDOCOUNT_SCALING=total,NULL_REPLICATES=${NULL_REPLICATES},RETAIN_UNCERTAINTY=0,UNCERTAINTY_SCALE=0" \
      "${REPO_ROOT}/analyses/microglia_less/run_paired_path_test.sbatch")
    merge_job=$(sbatch \
      --parsable \
      -p "${PARTITION}" \
      --time=00:40:00 \
      --dependency="afterok:${fit_job}" \
      --output="${LOG_ROOT}/full_testing_eb_${family}_a${strength}_merge.%j.out" \
      --error="${LOG_ROOT}/full_testing_eb_${family}_a${strength}_merge.%j.err" \
      --wrap="bash -lc 'module load ${PYTHON_MODULE} && source ${HOME}/venv/cpa_blackwell/bin/activate && PYTHONPATH=${REPO_ROOT} python ${REPO_ROOT}/extra_scripts/merge_paired_path_test.py --shards ${output_root}/shard_* --output-dir ${merged} --calibration empirical ${moderation_flag}'")
    merge_dependencies+=("${merge_job}")
    printf '%s\t%s\t%s\t%s\n' "${family}" "${strength}" "${fit_job}" "${merge_job}"
  done
  local dependency
  dependency=$(IFS=:; printf '%s' "${merge_dependencies[*]}")
  local select_job
  select_job=$(sbatch \
    --parsable \
    -p "${PARTITION}" \
    --dependency="afterok:${dependency}" \
    --output="${LOG_ROOT}/full_testing_eb_${family}_select.%j.out" \
    --error="${LOG_ROOT}/full_testing_eb_${family}_select.%j.err" \
    --export="ALL,CONFIG=${REPO_ROOT}/analyses/microglia_less/config.env,RESULT_ROOT=${result_root},STRENGTHS=${STRENGTHS}" \
    "${REPO_ROOT}/analyses/microglia_less/run_local_path_full_testing_eb_selection.sbatch")
  printf '%s\tselection\t%s\t%s\n' "${family}" "${select_job}" "${result_root}/selection"
}

printf 'family\tstrength\tjob\toutput\n'
if [[ "${RUN_FULL_OMNIBUS}" == 1 ]]; then
  submit_family omnibus_production "${RUN_ROOT}/differential/ec_block_full_omnibus_min4_candidates.pkl" ""
fi
if [[ "${RUN_FULL_CONDITION}" == 1 ]]; then
  submit_family condition "${RUN_ROOT}/differential/ec_block_condition_direct10_candidates.pkl" ""
fi
if [[ "${RUN_FULL_PAIRWISE}" == 1 ]]; then
  submit_family pairwise "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" "--moderate-variances"
fi
