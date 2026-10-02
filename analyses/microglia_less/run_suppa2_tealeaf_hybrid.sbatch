#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-191%48
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=suppa_tealeaf
#SBATCH --output=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/hybrid-%A_%a.out
#SBATCH --error=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/hybrid-%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

SHARD_COUNT=96
FOLD=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
OUT="${REPO_ROOT}/.cache/suppa2_statistics/tealeaf_hybrid_gene_wide/fold${FOLD}/shard_${SHARD}"

python "${REPO_ROOT}/extra_scripts/run_suppa2_tealeaf_hybrid.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold${FOLD}_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt" \
  --output-dir "${OUT}" \
  --shard-index "${SHARD}" \
  --shard-count "${SHARD_COUNT}" \
  --path-pseudocount 32 \
  --path-pseudocount-scaling total \
  --profile-event-mass --report-pseudocount 1 \
  --null-replicates 32
