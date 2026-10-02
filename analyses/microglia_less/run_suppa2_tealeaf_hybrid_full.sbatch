#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95%48
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=suppa_tealeaf_full
#SBATCH --output=.cache/suppa2_statistics/hybrid-full-%A_%a.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-full-%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

OUTPUT_TAG=${OUTPUT_TAG:-tealeaf_hybrid_full}
OUT="${REPO_ROOT}/.cache/suppa2_statistics/${OUTPUT_TAG}/shards/shard_${SLURM_ARRAY_TASK_ID}"
EXTRA_ARGS=()
if [[ "${SCAN_ALL_EVENTS:-1}" == "1" ]]; then
  EXTRA_ARGS+=(--scan-all-events)
fi
if [[ -n "${MAX_TESTS:-}" ]]; then
  EXTRA_ARGS+=(--max-tests "${MAX_TESTS}")
fi
python "${REPO_ROOT}/extra_scripts/run_suppa2_tealeaf_hybrid.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt" \
  --output-dir "${OUT}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --shard-count 96 \
  --path-pseudocount 64 \
  --path-pseudocount-scaling total \
  --profile-event-mass --report-pseudocount 1 \
  --null-replicates 32 \
  "${EXTRA_ARGS[@]}"
