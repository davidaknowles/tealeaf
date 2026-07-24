#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${EMBEDDING:?set EMBEDDING}"
: "${CELL_MAP:?set CELL_MAP}"
: "${LABELS:?set LABELS}"
: "${ELIGIBLE_CELLS:?set ELIGIBLE_CELLS}"
: "${SCORE_OUTPUT:?set SCORE_OUTPUT}"

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

group_args=()
if [[ -n "${GROUPS:-}" ]]; then
  group_args+=(--groups "${GROUPS}")
fi

python "${REPO_ROOT}/extra_scripts/score_reference_embedding.py" \
  --embedding "${EMBEDDING}" \
  --cell-map "${CELL_MAP}" \
  --labels "${LABELS}" \
  "${group_args[@]}" \
  --eligible-cells "${ELIGIBLE_CELLS}" \
  --name "${REFERENCE_NAME:-published_reference}" \
  --representation "${REFERENCE_REPRESENTATION:-published_embedding}" \
  --output-dir "${SCORE_OUTPUT}/${REFERENCE_NAME:-published_reference}"
