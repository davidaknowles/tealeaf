#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-63%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --job-name=event_input_controls
#SBATCH --output=.cache/event_information_20261007/fit_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/fit_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SHARD_COUNT=32
SOURCE_INDEX=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
SOURCES=(parsimony_binary original_binary)
SOURCE_NAME=${SOURCES[$SOURCE_INDEX]}
CONTROL="${REPO_ROOT}/.cache/event_information_20261007/${SOURCE_NAME}_paired"
OUTPUT="${REPO_ROOT}/.cache/event_information_20261007/full/${SOURCE_NAME}/shard_${SHARD_INDEX}"

# These are complete input controls, not a significance- or LR-ranked subset.
# Preserve the established full-data hybrid's concentration and reporting rule.
python "${REPO_ROOT}/extra_scripts/run_suppa2_tealeaf_hybrid.py" \
  --data-cache "${CONTROL}/prepared.pkl" \
  --features "${CONTROL}/features.txt" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --output-dir "${OUTPUT}" \
  --scan-all-events --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}" \
  --path-pseudocount 64 --path-pseudocount-scaling total \
  --profile-event-mass --report-pseudocount 1 --null-replicates 32
