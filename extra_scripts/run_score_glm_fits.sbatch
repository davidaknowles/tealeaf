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
fit_names=()
fit_prefixes=()
while IFS=$'\t' read -r name prefix; do
  [[ -n "${name}" && -n "${prefix}" ]] || continue
  fit_names+=("${name}")
  fit_prefixes+=("${prefix}")
done < "${FITS_FILE}"
if (( ${#fit_names[@]} == 0 )); then
  echo "No fits listed in ${FITS_FILE}" >&2
  exit 1
fi

fit_index=${FIT_INDEX:-${SLURM_ARRAY_TASK_ID:-}}
if [[ -n "${fit_index}" ]]; then
  if [[ ! "${fit_index}" =~ ^[0-9]+$ ]]; then
    echo "FIT_INDEX must be a nonnegative integer" >&2
    exit 1
  fi
  if (( fit_index < 0 || fit_index >= ${#fit_names[@]} )); then
    echo "FIT_INDEX ${fit_index} is outside the fit manifest" >&2
    exit 1
  fi
  fit_args+=(--fit "${fit_names[fit_index]}=${fit_prefixes[fit_index]}")
  output_dir="${SCORE_OUTPUT}/${fit_names[fit_index]}"
else
  for index in "${!fit_names[@]}"; do
    fit_args+=(--fit "${fit_names[index]}=${fit_prefixes[index]}")
  done
  output_dir="${SCORE_OUTPUT}"
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
  --output-dir "${output_dir}" \
  --device cuda
