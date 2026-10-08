#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --job-name=collate_sequence_null
#SBATCH --output=.cache/event_information_20261007/sequence_null_summary_%j.out

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
NULL_RUN=${SEQUENCE_NULL_RUN:-count_null}
RESULT=${SEQUENCE_NULL_RESULT:-event_sequence_count_null}
SOURCES=(parsimony_binary original_binary)
if [[ "${SEQUENCE_NULL_ORIGINAL_ONLY:-0}" == 1 ]]; then
  SOURCES=(original_binary)
fi
for SOURCE in "${SOURCES[@]}"; do
  for SCENARIO in strict biological nuisance; do
    python "${REPO_ROOT}/extra_scripts/summarize_ec_count_null.py" \
      --cache "${REPO_ROOT}/.cache/event_sequence_controls_L150_E001/${NULL_RUN}/${SOURCE}/${SCENARIO}" \
      --output-dir "${REPO_ROOT}/analyses/${RESULT}/${SOURCE}/${SCENARIO}"
  done
done
