#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00

set -euo pipefail
: "${JUNCTION_LIST:?set JUNCTION_LIST}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${LEAFCUTTER_ROOT:?set LEAFCUTTER_ROOT}"

module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
mkdir -p "${OUTPUT_DIR}"
python "${LEAFCUTTER_ROOT}/clustering/leafcutter_cluster_regtools.py" \
  --juncfiles "${JUNCTION_LIST}" --rundir "${OUTPUT_DIR}" \
  --outprefix benchmark --maxintronlen 500000 --minclureads 30 \
  --mincluratio 0.001
test -s "${OUTPUT_DIR}/benchmark_perind_numers.counts.gz"
