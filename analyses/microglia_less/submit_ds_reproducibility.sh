#!/bin/bash
set -euo pipefail

source "$(dirname "$0")/config.env"
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

BENCHMARK_ROOT="${RUN_ROOT}/junction_benchmark"
BUNDLE="${BENCHMARK_ROOT}/pseudobulk_junctions"
REPRO_ROOT="${BENCHMARK_ROOT}/reproducibility"
SUBJECT_FOLDS="${REPRO_ROOT}/subject_folds.tsv"
LEAFCUTTER_ROOT=${LEAFCUTTER_ROOT:-${HOME}/projects/leafcutter}
SCQUINT_ROOT=${SCQUINT_ROOT:-${HOME}/projects/scquint}

PYTHON="${HOME}/venv/cpa_blackwell/bin/python"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${REPRO_ROOT}/logs"

"${PYTHON}" "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" split-subjects \
  --bundle "${BUNDLE}" --folds 2 --seed 20260805 --output "${SUBJECT_FOLDS}"

job_file="${REPRO_ROOT}/jobs.tsv"
printf 'fold\tmethod\tjob_id\n' > "${job_file}"
for fold in 0 1; do
  fold_root="${REPRO_ROOT}/fold${fold}"
  contrasts="${fold_root}/contrasts.json"
  mkdir -p "${fold_root}/logs"
  "${PYTHON}" "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" plan \
    --bundle "${BUNDLE}" --subject-folds "${SUBJECT_FOLDS}" \
    --subject-fold "${fold}" --effect cell_type --min-subjects 4 \
    --output "${contrasts}"
  n=$("${PYTHON}" -c 'import json,sys; print(len(json.load(open(sys.argv[1]))))' "${contrasts}")

  leaf_job=$(sbatch --parsable --array="0-$((n - 1))%32" \
    --job-name="dsrep_leaf_f${fold}" \
    --output="${fold_root}/logs/leaf.%A_%a.out" \
    --export="ALL,REPO_ROOT=${REPO_ROOT},BUNDLE=${BUNDLE},CONTRASTS=${contrasts},COUNTS=${BENCHMARK_ROOT}/leafcutter/clustering/benchmark_perind_numers.counts.gz,OUTPUT_DIR=${fold_root}/leafcutter/tests,LEAFCUTTER_ROOT=${LEAFCUTTER_ROOT},PYTHON_MODULE=${PYTHON_MODULE}" \
    "${REPO_ROOT}/extra_scripts/run_leafcutter_contrast.sbatch")

  scquint_job=$(sbatch --parsable --array="0-$((n - 1))%32" \
    --job-name="dsrep_scq_f${fold}" \
    --output="${fold_root}/logs/scquint.%A_%a.out" \
    --export="ALL,REPO_ROOT=${REPO_ROOT},BUNDLE=${BUNDLE},CONTRASTS=${contrasts},OUTPUT_DIR=${fold_root}/scquint/tests,SCQUINT_ROOT=${SCQUINT_ROOT},PYTHON_MODULE=${PYTHON_MODULE}" \
    "${REPO_ROOT}/extra_scripts/run_scquint_contrast.sbatch")

  majiq_job=$(sbatch --parsable --array="0-$((n - 1))%16" \
    --job-name="dsrep_majiq_f${fold}" \
    --output="${fold_root}/logs/majiq.%A_%a.out" \
    --export="ALL,REPO_ROOT=${REPO_ROOT},CONTRASTS=${contrasts},PSICOV=${BENCHMARK_ROOT}/majiq/all.psicov,SPLICEGRAPH=${BENCHMARK_ROOT}/majiq/combined.splicegraph,ANNOTATED_SPLICEGRAPH=${BENCHMARK_ROOT}/majiq/annotated.splicegraph,OUTPUT_DIR=${fold_root}/majiq/tests,MAJIQ_LICENSE_FILE=${MAJIQ_LICENSE_FILE}" \
    "${REPO_ROOT}/extra_scripts/run_majiq_contrast.sbatch")

  summary_job=$(sbatch --parsable \
    --dependency="afterok:${leaf_job}:${scquint_job}:${majiq_job}" \
    --job-name="dsrep_summary_f${fold}" \
    --output="${fold_root}/logs/summary.%j.out" \
    --export="ALL,REPO_ROOT=${REPO_ROOT},BENCHMARK_ROOT=${fold_root},PYTHON_MODULE=${PYTHON_MODULE}" \
    "${REPO_ROOT}/extra_scripts/run_junction_summary.sbatch")

  output_tag="ec_block_reproducibility_fold${fold}"
  candidate_cache="${RUN_ROOT}/differential/${output_tag}_candidates.pkl"
  if [[ ! -s "${candidate_cache}" ]]; then
    "${PYTHON}" "${REPO_ROOT}/extra_scripts/run_ec_block_glmm.py" \
      --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
      --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt" \
      --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz" \
      --output-dir "${REPRO_ROOT}/candidate_screen_fold${fold}" \
      --methods laplace_multinomial --calibration lrt_bic \
      --candidate-cache "${candidate_cache}" --test-effect cell_type \
      --min-gene-umis 25 --min-celltype-mice 4 \
      --subject-folds "${SUBJECT_FOLDS}" --subject-fold "${fold}" \
      --shard-index 0 --shard-count 1 --dry-run
  fi
  tealeaf_job=$(sbatch --parsable --array="0-95%48" \
    --job-name="dsrep_tealeaf_f${fold}" \
    --export="ALL,CONFIG=${REPO_ROOT}/analyses/microglia_less/config.env,OUTPUT_TAG=${output_tag},CANDIDATE_CACHE=${candidate_cache},CALIBRATION=lrt_bic,METHODS=laplace_multinomial,NULL_REPLICATES=0,TEST_EFFECT=cell_type,MIN_GENE_UMIS=25,MIN_CELLTYPE_MICE=4,SUBJECT_FOLDS=${SUBJECT_FOLDS},SUBJECT_FOLD=${fold}" \
    "${REPO_ROOT}/analyses/microglia_less/run_ec_block_glmm.sbatch")

  tealeaf_output="${fold_root}/tealeaf"
  merge_command="source ${HOME}/venv/cpa_blackwell/bin/activate && PYTHONPATH=${REPO_ROOT} python ${REPO_ROOT}/extra_scripts/merge_ec_block_asymptotic.py --shards ${RUN_ROOT}/differential/${output_tag}_shard* --candidate-cache ${candidate_cache} --output-dir ${tealeaf_output}"
  merge_job=$(sbatch --parsable --dependency="afterok:${tealeaf_job}" \
    --job-name="dsrep_merge_f${fold}" \
    --output="${fold_root}/logs/tealeaf_merge.%j.out" \
    --wrap="bash -lc '${merge_command}'")

  printf '%s\t%s\t%s\n' "${fold}" leafcutter "${leaf_job}" >> "${job_file}"
  printf '%s\t%s\t%s\n' "${fold}" scquint "${scquint_job}" >> "${job_file}"
  printf '%s\t%s\t%s\n' "${fold}" majiq "${majiq_job}" >> "${job_file}"
  printf '%s\t%s\t%s\n' "${fold}" external_summary "${summary_job}" >> "${job_file}"
  printf '%s\t%s\t%s\n' "${fold}" tealeaf "${tealeaf_job}" >> "${job_file}"
  printf '%s\t%s\t%s\n' "${fold}" tealeaf_merge "${merge_job}" >> "${job_file}"
done

cat "${job_file}"
