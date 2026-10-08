#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --job-name=sequence_event_null
#SBATCH --output=.cache/event_information_20261007/sequence_null_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/sequence_null_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SHARD_COUNT=16
EVENTS=${SEQUENCE_NULL_EVENTS:-32}
DRAWS=${SEQUENCE_NULL_DRAWS:-2}
NULL_RUN=${SEQUENCE_NULL_RUN:-count_null}
SOURCE_INDEX=$((SLURM_ARRAY_TASK_ID / (3 * SHARD_COUNT)))
SCENARIO_INDEX=$(((SLURM_ARRAY_TASK_ID / SHARD_COUNT) % 3))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
SOURCES=(parsimony_binary original_binary)
SCENARIOS=(strict biological nuisance)
SOURCE=${SOURCES[$SOURCE_INDEX]}
SCENARIO=${SCENARIOS[$SCENARIO_INDEX]}
EXTRA=()
if [[ "${SCENARIO}" == biological ]]; then
  EXTRA+=(--residual-concentration 20)
elif [[ "${SCENARIO}" == nuisance ]]; then
  EXTRA+=(--within-path-type-scale 1)
fi
python "${REPO_ROOT}/extra_scripts/audit_hybrid_input_count_null.py" \
  --cache "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001" \
  --family-cache "${REPO_ROOT}/.cache/event_information_20261007" \
  --source "${SOURCE}" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold0_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --output-dir "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001/${NULL_RUN}/${SOURCE}/${SCENARIO}/shard_${SHARD_INDEX}" \
  --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}" \
  --events "${EVENTS}" --draws "${DRAWS}" --subject-scale .5 --inference mixed-score \
  --max-iter 2000 --null-multistart --information-metric reference --scalar-fast \
  "${EXTRA[@]}"
