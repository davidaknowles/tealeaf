#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-63%32
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=hybrid_split_fix
#SBATCH --output=.cache/suppa2_statistics/hybrid-split-fix-%A_%a.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-split-fix-%A_%a.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SHARD_COUNT=32
FOLD=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
DEBUG_ROOT=${DEBUG_ROOT:-.cache/suppa2_statistics/hybrid_lr_analytic}
python extra_scripts/run_suppa2_tealeaf_hybrid.py \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold${FOLD}_candidates.pkl" \
  --event-catalog analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/full_data/event_catalog.tsv.gz \
  --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt" \
  --output-dir "${DEBUG_ROOT}/profile32_split/fold${FOLD}/shard_${SHARD}" \
  --shard-index "$SHARD" --shard-count "$SHARD_COUNT" \
  --path-pseudocount 32 --path-pseudocount-scaling total \
  --profile-event-mass --report-pseudocount 1 --null-replicates 32
