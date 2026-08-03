#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=08:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${RUN_LIST:?tab-separated run and STAR output directory}"
: "${CELL_METADATA:?set CELL_METADATA}"

line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${RUN_LIST}")
IFS=$'\t' read -r run output_dir <<< "${line}"
if [[ -z "${run}" || -z "${output_dir}" ]]; then
  echo "No run at array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python "${REPO_ROOT}/extra_scripts/tag_spliced_pseudobulk_bam.py" \
  --input "${output_dir}/Aligned.sortedByCoord.out.bam" \
  --cell-metadata "${CELL_METADATA}" --run "${run}" \
  --output "${output_dir}/spliced_pseudobulk.bam"
