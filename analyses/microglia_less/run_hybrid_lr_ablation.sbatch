#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-63%32
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=hybrid_lr_ablation
#SBATCH --output=.cache/suppa2_statistics/hybrid-lr-ablation-%A_%a.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-lr-ablation-%A_%a.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SHARD_COUNT=32
SHARD=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
VARIANT=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
OPTIONS=()
TEST_ZETA=${TEST_ZETA:-64}
DEBUG_ROOT=${DEBUG_ROOT:-.cache/suppa2_statistics/hybrid_lr_analytic}
REPORT_ZETA=${REPORT_ZETA:-1}
REPORT_ARGS=(--report-pseudocount "$REPORT_ZETA")
TAG=fixed${TEST_ZETA}_report1
if [[ "$REPORT_ZETA" == none ]]; then
  REPORT_ARGS=()
  TAG=fixed${TEST_ZETA}_testonly
fi
if [[ "$VARIANT" == 1 ]]; then
  TAG=profile${TEST_ZETA}_report1
  OPTIONS+=(--profile-event-mass)
fi
python extra_scripts/run_suppa2_tealeaf_hybrid.py \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --event-catalog analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/full_data/event_catalog.tsv.gz \
  --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt" \
  --output-dir "${DEBUG_ROOT}/${TAG}/shards/shard_${SHARD}" \
  --shard-index "$SHARD" --shard-count "$SHARD_COUNT" \
  --path-pseudocount "$TEST_ZETA" --path-pseudocount-scaling total \
  "${REPORT_ARGS[@]}" --export-path-usage --null-replicates 32 \
  --scan-all-events "${OPTIONS[@]}"
