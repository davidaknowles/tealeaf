#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-15%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --job-name=sequence_observed_score
#SBATCH --output=.cache/event_information_20261007/sequence_observed_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/sequence_observed_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SHARD_COUNT=8
SOURCE_INDEX=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
SOURCES=(parsimony_binary original_binary)
SOURCE=${SOURCES[$SOURCE_INDEX]}
python "${REPO_ROOT}/extra_scripts/audit_observed_event_score.py" \
  --cache "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001" --source "${SOURCE}" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --native-support "${REPO_ROOT}/analyses/event_information_audit/native_event_support.tsv.gz" \
  --null-panel "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001/count_null/${SOURCE}/strict/shard_0/settings.json" \
  --output-dir "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001/observed/${SOURCE}/shard_${SHARD_INDEX}" \
  --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}" --information-metric reference
