#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${BUNDLE:?set BUNDLE}"
: "${CONTRASTS:?set CONTRASTS}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${SCQUINT_ROOT:?set SCQUINT_ROOT}"
index=${SLURM_ARRAY_TASK_ID}

module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${OUTPUT_DIR}"
python "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" scquint \
  --bundle "${BUNDLE}" --contrasts "${CONTRASTS}" --contrast-index "${index}" \
  --scquint-root "${SCQUINT_ROOT}" --jobs "${SLURM_CPUS_PER_TASK:-4}" \
  --output "${OUTPUT_DIR}/${index}.tsv"
