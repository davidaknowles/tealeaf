#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=assess_event_mass_null
#SBATCH --output=.cache/event_information_20261007/assess_event_mass_%j.out

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

for SCENARIO in mass_only mass_and_within; do
  python "${REPO_ROOT}/extra_scripts/summarize_ec_count_null.py" \
    --cache "${REPO_ROOT}/.cache/event_sequence_event_mass_null/original_binary/${SCENARIO}" \
    --output-dir "${REPO_ROOT}/analyses/event_sequence_event_mass_null/original_binary/${SCENARIO}"
done
