#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=4-00:00:00

set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${ALEVIN_DIR:?set ALEVIN_DIR}"
: "${SALMON_REF:?set SALMON_REF}"
: "${RUN_DIR:?set RUN_DIR}"
: "${DESIGN:?set DESIGN to binary or weighted}"
: "${METHOD:?set METHOD}"
module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${RUN_DIR}"

CV_CELLS=${CV_CELLS:-0}
MIN_CELL_UMIS=${MIN_CELL_UMIS:-500}
MIN_HALF_UMIS=${MIN_HALF_UMIS:-500}
PRIMER_SAMPLING_MODEL=${PRIMER_SAMPLING_MODEL:-oligodt_tpm}
BATCH_CELLS=${BATCH_CELLS:-16384}
MAX_RANK=${MAX_RANK:-512}
EXACT_INNER_STEPS=${EXACT_INNER_STEPS:-32}
ANALYSIS_TAG=${ANALYSIS_TAG:-standard}
REPORT="${RUN_DIR}/cv_${ANALYSIS_TAG}_${DESIGN}_${METHOD}_theta.json"

multiplier_args=()
rank_args=()
if [[ "${METHOD}" == factorized ]]; then
  for value in 1 2 4 8 16 32 64 128; do
    rank_args+=(--rank-candidate "${value}")
  done
  max_iter=${MAX_ITER:-4096}
  expansion_factor=2
  selection_rule=one_standard_error
  cv_tol=1e-5
elif [[ "${METHOD}" == admm_factorized ]]; then
  for value in 0 1e-9 1e-8 1e-7 1e-6 1e-5; do
    multiplier_args+=(--multiplier "${value}")
  done
  max_iter=${MAX_ITER:-4096}
  expansion_factor=10
  selection_rule=minimum
  cv_tol=1e-5
elif [[ "${METHOD}" == frank_wolfe_penalized ]]; then
  for value in 0 0.25 1 4 16 64 256; do
    multiplier_args+=(--multiplier "${value}")
  done
  max_iter=${MAX_ITER:-2048}
  expansion_factor=4
  selection_rule=one_se_variance_retention
  cv_tol=1e-4
else
  echo "Unsupported scalable method: ${METHOD}" >&2
  exit 1
fi

primer_args=()
if [[ -n "${PRIMER_PAIRS:-}" ]]; then
  primer_args+=(
    --primer-pairs "${PRIMER_PAIRS}"
    --min-half-umis "${MIN_HALF_UMIS}"
    --primer-sampling-model "${PRIMER_SAMPLING_MODEL}"
  )
fi

python "${REPO_ROOT}/extra_scripts/tune_single_cell_glm.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}" \
  --design "${DESIGN}" \
  --method "${METHOD}" \
  --output "${REPORT}" \
  --cells "${CV_CELLS}" \
  --min-cell-umis "${MIN_CELL_UMIS}" \
  --folds 3 \
  --max-iter "${max_iter}" \
  --rank 64 \
  --max-rank "${MAX_RANK}" \
  --power-iter 5 \
  --scale-power-iter 10 \
  --batch-cells "${BATCH_CELLS}" \
  --data-backend auto \
  --polish-max-iter 32 \
  --exact-inner-steps "${EXACT_INNER_STEPS}" \
  --device cuda \
  --max-grid-expansions 4 \
  --grid-expansion-factor "${expansion_factor}" \
  --selection-rule "${selection_rule}" \
  --require-converged \
  --require-nondegenerate \
  --tol "${cv_tol}" \
  --profile-variance-retention 0.9 \
  "${primer_args[@]}" \
  "${multiplier_args[@]}" \
  "${rank_args[@]}"
