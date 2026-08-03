#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${BUNDLE:?set BUNDLE}"
: "${CONTRASTS:?set CONTRASTS}"
: "${COUNTS:?set COUNTS}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${LEAFCUTTER_ROOT:?set LEAFCUTTER_ROOT}"
index=${SLURM_ARRAY_TASK_ID}
mkdir -p "${OUTPUT_DIR}/metadata" "${OUTPUT_DIR}/raw" "${OUTPUT_DIR}/normalized"

module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" leafcutter-metadata \
  --bundle "${BUNDLE}" --contrasts "${CONTRASTS}" --contrast-index "${index}" \
  --output "${OUTPUT_DIR}/metadata/${index}.tsv"

module purge
module load leafcutter/0.2.9
Rscript "${LEAFCUTTER_ROOT}/scripts/leafcutter_ds.R" \
  --num_threads "${SLURM_CPUS_PER_TASK:-4}" \
  --min_samples_per_intron 3 --min_samples_per_group 3 --min_coverage 20 \
  --output_prefix "${OUTPUT_DIR}/raw/${index}" \
  "${COUNTS}" "${OUTPUT_DIR}/metadata/${index}.tsv"

module purge
module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" normalize-leafcutter \
  --contrasts "${CONTRASTS}" --contrast-index "${index}" \
  --input "${OUTPUT_DIR}/raw/${index}_cluster_significance.txt" \
  --output "${OUTPUT_DIR}/normalized/${index}.tsv"
