#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-191%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=event_score_reference
#SBATCH --output=.cache/event_information_20261007/reference_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/reference_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SHARD_COUNT=64
COHORTS=(full fold0 fold1)
COHORT=${COHORTS[$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))]}
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
python "${REPO_ROOT}/extra_scripts/reassess_event_score_archive.py" \
  --input-dir "${REPO_ROOT}/.cache/event_information_20261007/mixed_score/original_binary/${COHORT}/shard_${SHARD_INDEX}" \
  --output-dir "${REPO_ROOT}/.cache/event_information_20261007/mixed_score_reference/original_binary/${COHORT}/shard_${SHARD_INDEX}" \
  --information-metric reference --refit-missing
