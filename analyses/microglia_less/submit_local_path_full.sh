#!/bin/bash
set -euo pipefail

source "$(dirname "$0")/config.env"

PARTITION=${PARTITION:-cpu}
TIME_LIMIT=${TIME_LIMIT:-02:00:00}
SHARD_LIMIT=${SHARD_LIMIT:-48}
NULL_REPLICATES=${NULL_REPLICATES:-32}
LOG_ROOT="${RUN_ROOT}/junction_benchmark/reproducibility/logs"

submit_analysis() {
  local contrast=$1
  local candidates=$2
  local uncertainty_scale=$3
  local moderation_flag=$4
  local output_root="${RUN_ROOT}/differential/local_path_full_${contrast}_production_shards"
  local merged="${RUN_ROOT}/differential/local_path_full_${contrast}_production"
  local fit_job
  local merge_job
  fit_job=$(sbatch --parsable -p "${PARTITION}" --time="${TIME_LIMIT}" --array="0-95%${SHARD_LIMIT}" --output="${LOG_ROOT}/full_${contrast}_production.%A_%a.out" --error="${LOG_ROOT}/full_${contrast}_production.%A_%a.err" --export="ALL,CONFIG=${REPO_ROOT}/analyses/microglia_less/config.env,CANDIDATE_CACHE=${candidates},OUTPUT_ROOT=${output_root},PATH_PSEUDOCOUNT=2,PATH_PRIOR_CENTER=uniform,NULL_REPLICATES=${NULL_REPLICATES},RETAIN_UNCERTAINTY=0,UNCERTAINTY_SCALE=${uncertainty_scale}" "${REPO_ROOT}/analyses/microglia_less/run_paired_path_test.sbatch")
  merge_job=$(sbatch --parsable -p "${PARTITION}" --time=00:40:00 --dependency="afterok:${fit_job}" --output="${LOG_ROOT}/full_${contrast}_production_merge.%j.out" --error="${LOG_ROOT}/full_${contrast}_production_merge.%j.err" --wrap="bash -lc 'module load ${PYTHON_MODULE} && source ${HOME}/venv/cpa_blackwell/bin/activate && PYTHONPATH=${REPO_ROOT} python ${REPO_ROOT}/extra_scripts/merge_paired_path_test.py --shards ${output_root}/shard_* --output-dir ${merged} ${moderation_flag} --calibration empirical'")
  printf '%s\t%s\t%s\n' "${contrast}" "${fit_job}" "${merge_job}"
}

printf 'contrast\tfit_job\tmerge_job\n'
submit_analysis celltype "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" 1 "--moderate-variances"
submit_analysis condition "${RUN_ROOT}/differential/ec_block_condition_direct10_candidates.pkl" 0 ""
