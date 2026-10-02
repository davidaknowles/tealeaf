#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --job-name=suppa_hybrid_reference
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/hybrid_reference_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/hybrid_reference_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REPRO="${RUN_ROOT}/junction_benchmark/reproducibility"
HYBRID="${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/split_data"
python "${REPO_ROOT}/extra_scripts/evaluate_suppa2_statistics.py" \
  --tealeaf-tests "${REPRO}/fold0/tealeaf_paired_path_total_a32_production/paired_path.tsv" \
  --tealeaf-tests "${REPRO}/fold1/tealeaf_paired_path_total_a32_production/paired_path.tsv" \
  --comparison "SUPPA2 definitions; Tealeaf EC; current reference" "${HYBRID}/fold0_tests.tsv.gz" "${HYBRID}/fold1_tests.tsv.gz" \
  --output-dir "${HYBRID}/production_reference"
