#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-47%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=event_mixed_count_null
#SBATCH --output=.cache/event_information_20261007/mixed_null_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/mixed_null_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SHARD_COUNT=16
SCENARIO_INDEX=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
SCENARIOS=(strict biological nuisance)
SCENARIO=${SCENARIOS[$SCENARIO_INDEX]}
SCORE_COORDINATE=${SCORE_COORDINATE:-ilr}
SCORE_SOURCE=${SCORE_SOURCE:-original_binary}
SCORE_INFORMATION_METRIC=${SCORE_INFORMATION_METRIC:-absolute}
SCORE_EVENTS=${SCORE_EVENTS:-32}
SCORE_DRAWS=${SCORE_DRAWS:-2}
SCORE_COUNT_LIKELIHOOD=${SCORE_COUNT_LIKELIHOOD:-multinomial}
SCORE_EC_OPPORTUNITY_SCALE=${SCORE_EC_OPPORTUNITY_SCALE:-0}
SCORE_KERNEL_UNITS=${SCORE_KERNEL_UNITS:-prepared}
SCORE_SIMULATION_KERNEL_UNITS=${SCORE_SIMULATION_KERNEL_UNITS:-analysis}
if [[ "${SCORE_INFORMATION_METRIC}" != "absolute" ]]; then
  SCORE_OUTPUT="${REPO_ROOT}/.cache/event_information_20261007/count_null_score_${SCORE_COORDINATE}_${SCORE_INFORMATION_METRIC}/${SCORE_SOURCE}"
elif [[ "${SCORE_COORDINATE}" == "ilr" && "${SCORE_SOURCE}" == "original_binary" ]]; then
  SCORE_OUTPUT="${REPO_ROOT}/.cache/event_information_20261007/count_null_mixed_original"
else
  SCORE_OUTPUT="${REPO_ROOT}/.cache/event_information_20261007/count_null_score_${SCORE_COORDINATE}/${SCORE_SOURCE}"
fi
if [[ "${SCORE_EVENTS}" != "32" || "${SCORE_DRAWS}" != "2" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_events${SCORE_EVENTS}_draws${SCORE_DRAWS}"
fi
if [[ "${SCORE_COUNT_LIKELIHOOD}" != "multinomial" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_${SCORE_COUNT_LIKELIHOOD}"
fi
if [[ "${SCORE_EC_OPPORTUNITY_SCALE}" != "0" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_opportunity${SCORE_EC_OPPORTUNITY_SCALE}"
fi
if [[ "${SCORE_KERNEL_UNITS}" != "prepared" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_analysis${SCORE_KERNEL_UNITS}"
fi
if [[ "${SCORE_SIMULATION_KERNEL_UNITS}" != "analysis" ]]; then
  SCORE_OUTPUT="${SCORE_OUTPUT}_simulation${SCORE_SIMULATION_KERNEL_UNITS}"
fi
EXTRA=()
if [[ "${SCENARIO}" == "biological" ]]; then
  EXTRA+=(--residual-concentration 20)
elif [[ "${SCENARIO}" == "nuisance" ]]; then
  EXTRA+=(--within-path-type-scale 1)
fi
python "${REPO_ROOT}/extra_scripts/audit_hybrid_input_count_null.py" \
  --cache "${REPO_ROOT}/.cache/event_information_20261007" --source "${SCORE_SOURCE}" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold0_candidates.pkl" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --output-dir "${SCORE_OUTPUT}/${SCENARIO}/shard_${SHARD_INDEX}" \
  --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}" --events "${SCORE_EVENTS}" --draws "${SCORE_DRAWS}" \
  --inference mixed-score --score-coordinate "${SCORE_COORDINATE}" --information-metric "${SCORE_INFORMATION_METRIC}" --count-likelihood "${SCORE_COUNT_LIKELIHOOD}" --ec-opportunity-scale "${SCORE_EC_OPPORTUNITY_SCALE}" --kernel-units "${SCORE_KERNEL_UNITS}" --simulation-kernel-units "${SCORE_SIMULATION_KERNEL_UNITS}" --max-iter 2000 --null-multistart "${EXTRA[@]}"
