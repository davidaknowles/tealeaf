#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-15%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=split_directions
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/split_coverage_audit/directions_%A_%a.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/split_coverage_audit/directions_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
fold=$((SLURM_ARRAY_TASK_ID / 8))
shard=$((SLURM_ARRAY_TASK_ID % 8))
python "${REPO_ROOT}/extra_scripts/recover_archived_path_directions.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold${fold}_candidates.pkl" \
  --output-dir "${REPO_ROOT}/.cache/split_coverage_audit/archived_directions/fold${fold}/shard_${shard}" \
  --null-replicates 0 --path-pseudocount 32 --path-prior-center uniform \
  --path-pseudocount-scaling total --shard-index "${shard}" --shard-count 8
