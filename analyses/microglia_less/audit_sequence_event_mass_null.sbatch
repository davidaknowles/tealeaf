#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-31%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --job-name=sequence_event_mass_null
#SBATCH --output=.cache/event_information_20261007/event_mass_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/event_mass_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SHARD_COUNT=16
SCENARIOS=(mass_only mass_and_within)
SCENARIO=${SCENARIOS[$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))]}
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
EXTRA=()
if [[ "${SCENARIO}" == mass_and_within ]]; then
  EXTRA+=(--within-path-type-scale 1)
fi

python "${REPO_ROOT}/extra_scripts/audit_hybrid_input_count_null.py" \
  --cache "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001" \
  --family-cache "${REPO_ROOT}/.cache/event_information_20261007" --source original_binary \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold0_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --output-dir "${REPO_ROOT}/.cache/event_sequence_event_mass_null/original_binary/${SCENARIO}/shard_${SHARD_INDEX}" \
  --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}" --events 128 --draws 2 \
  --inference mixed-score --max-iter 2000 --null-multistart --scalar-fast \
  --information-metric reference --event-mass-type-scale 2 "${EXTRA[@]}"
