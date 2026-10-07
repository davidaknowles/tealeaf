#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-15%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=observed_event_score
#SBATCH --output=.cache/event_information_20261007/observed_coordinate_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/observed_coordinate_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SCORE_COORDINATE=${SCORE_COORDINATE:-proportion}
SCORE_INFORMATION_METRIC=${SCORE_INFORMATION_METRIC:-absolute}
SCORE_COUNT_LIKELIHOOD=${SCORE_COUNT_LIKELIHOOD:-multinomial}
SCORE_KERNEL_UNITS=${SCORE_KERNEL_UNITS:-prepared}
SCORE_OUTPUT="${REPO_ROOT}/.cache/event_information_20261007/observed_score_${SCORE_COORDINATE}"
if [[ "${SCORE_INFORMATION_METRIC}" != "absolute" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_${SCORE_INFORMATION_METRIC}"
fi
if [[ "${SCORE_COUNT_LIKELIHOOD}" != "multinomial" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_${SCORE_COUNT_LIKELIHOOD}"
fi
if [[ "${SCORE_KERNEL_UNITS}" != "prepared" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_${SCORE_KERNEL_UNITS}"
fi
SHARD_COUNT=8
SOURCE_INDEX=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
SOURCES=(parsimony_binary original_binary)
SOURCE_NAME=${SOURCES[$SOURCE_INDEX]}
python "${REPO_ROOT}/extra_scripts/audit_observed_event_score.py" \
  --cache "${REPO_ROOT}/.cache/event_information_20261007" --source "${SOURCE_NAME}" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --native-support "${REPO_ROOT}/analyses/event_information_audit/native_event_support.tsv.gz" \
  --null-panel "${REPO_ROOT}/.cache/event_information_20261007/count_null/parsimony_binary/strict/shard_0/settings.json" \
  --output-dir "${SCORE_OUTPUT}/${SOURCE_NAME}/shard_${SHARD_INDEX}" \
  --score-coordinate "${SCORE_COORDINATE}" --information-metric "${SCORE_INFORMATION_METRIC}" --count-likelihood "${SCORE_COUNT_LIKELIHOOD}" --kernel-units "${SCORE_KERNEL_UNITS}" --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}"
