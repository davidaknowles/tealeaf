#!/bin/bash
set -euo pipefail

source "$(dirname "$0")/config.env"

BENCHMARK_ROOT="${RUN_ROOT}/junction_benchmark"
REPRO_ROOT="${BENCHMARK_ROOT}/reproducibility"
PARTITION=${PARTITION:-cpu}
TIME_LIMIT=${TIME_LIMIT:-02:00:00}
PATH_PSEUDOCOUNT=${PATH_PSEUDOCOUNT:-2}
OUTPUT_TAG=${OUTPUT_TAG:-tealeaf_paired_path}
job_file="${REPRO_ROOT}/paired_path_jobs.tsv"
mkdir -p "${REPRO_ROOT}/logs"
printf 'fold\tstage\tjob_id\n' > "${job_file}"
merge_jobs=()

for fold in 0 1; do
  fold_root="${REPRO_ROOT}/fold${fold}"
  mkdir -p "${fold_root}/logs"
  output_root="${fold_root}/${OUTPUT_TAG}_shards"
  candidate_cache="${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold${fold}_candidates.pkl"
  fit_job=$(sbatch --parsable -p "${PARTITION}" --time="${TIME_LIMIT}" --array=0-95%48 --output="${fold_root}/logs/${OUTPUT_TAG}.%A_%a.out" --error="${fold_root}/logs/${OUTPUT_TAG}.%A_%a.err" --export="ALL,CONFIG=${REPO_ROOT}/analyses/microglia_less/config.env,CANDIDATE_CACHE=${candidate_cache},OUTPUT_ROOT=${output_root},PATH_PSEUDOCOUNT=${PATH_PSEUDOCOUNT},NULL_REPLICATES=32" "${REPO_ROOT}/analyses/microglia_less/run_paired_path_test.sbatch")
  merge_command="source ${HOME}/venv/cpa_blackwell/bin/activate && PYTHONPATH=${REPO_ROOT} python ${REPO_ROOT}/extra_scripts/merge_paired_path_test.py --shards ${output_root}/shard_* --output-dir ${fold_root}/${OUTPUT_TAG} --moderate-variances --calibration empirical"
  merge_job=$(sbatch --parsable -p "${PARTITION}" --time="${TIME_LIMIT}" --dependency="afterok:${fit_job}" --output="${fold_root}/logs/${OUTPUT_TAG}_merge.%j.out" --error="${fold_root}/logs/${OUTPUT_TAG}_merge.%j.err" --wrap="bash -lc '${merge_command}'")
  merge_jobs+=("${merge_job}")
  printf '%s\t%s\t%s\n' "${fold}" fit "${fit_job}" >> "${job_file}"
  printf '%s\t%s\t%s\n' "${fold}" merge "${merge_job}" >> "${job_file}"
done

score_command="source ${HOME}/venv/cpa_blackwell/bin/activate && PYTHONPATH=${REPO_ROOT} python ${REPO_ROOT}/extra_scripts/score_ds_reproducibility.py --fold-dir ${REPRO_ROOT}/fold0 --fold-dir ${REPRO_ROOT}/fold1 --junction-table ${BENCHMARK_ROOT}/pseudobulk_junctions_junctions.tsv.gz --leafcutter-counts ${BENCHMARK_ROOT}/leafcutter/clustering/benchmark_perind_numers.counts.gz --leafcutter-map ${REPRO_ROOT}/leafcutter_cluster_gene.tsv.gz --tealeaf-subdir ${OUTPUT_TAG} --comparison-input-subdir comparison_majiq_min3_cov3 --tealeaf-label 'Tealeaf direct local path' --tealeaf-format paired_path --min-paired-subjects 4 --match-pairs --output-dir ${REPRO_ROOT}/${OUTPUT_TAG}_comparison"
score_job=$(sbatch --parsable -p "${PARTITION}" --time="${TIME_LIMIT}" --dependency="afterok:${merge_jobs[0]}:${merge_jobs[1]}" --output="${REPRO_ROOT}/logs/${OUTPUT_TAG}_score.%j.out" --error="${REPRO_ROOT}/logs/${OUTPUT_TAG}_score.%j.err" --wrap="bash -lc '${score_command}'")
printf '%s\t%s\t%s\n' all score "${score_job}" >> "${job_file}"
cat "${job_file}"
