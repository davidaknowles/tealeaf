#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00

set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${FITS_FILE:?set FITS_FILE to a tab-separated name/prefix manifest}"
: "${LABELS:?set LABELS}"
: "${TRANSCRIPT_TO_GENE:?set TRANSCRIPT_TO_GENE}"
: "${SCORE_OUTPUT:?set SCORE_OUTPUT}"

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

fit_args=()
while IFS=$'\t' read -r name prefix; do
  [[ -n "${name}" && -n "${prefix}" ]] || continue
  fit_args+=(--fit "${name}=${prefix}")
done < "${FITS_FILE}"
if (( ${#fit_args[@]} == 0 )); then
  echo "No fits listed in ${FITS_FILE}" >&2
  exit 1
fi

group_args=()
if [[ -n "${GROUPS:-}" ]]; then
  group_args+=(--groups "${GROUPS}")
fi

python "${REPO_ROOT}/extra_scripts/score_single_cell_glm_labels.py" \
  "${fit_args[@]}" \
  --labels "${LABELS}" \
  "${group_args[@]}" \
  --transcript-to-gene "${TRANSCRIPT_TO_GENE}" \
  --output-dir "${SCORE_OUTPUT}" \
  --device cuda
