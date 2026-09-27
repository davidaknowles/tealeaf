#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=merge_suppa_tealeaf
#SBATCH --output=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/hybrid-merge-%j.out
#SBATCH --error=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/hybrid-merge-%j.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

BASE="${REPO_ROOT}/.cache/suppa2_statistics/tealeaf_hybrid_gene_wide"
for fold in 0 1; do
  python "${REPO_ROOT}/extra_scripts/merge_paired_path_test.py" \
    --shards "${BASE}/fold${fold}"/shard_* \
    --output-dir "${BASE}/fold${fold}_merged" \
    --calibration empirical
done

RUN_REPRO="${RUN_ROOT}/junction_benchmark/reproducibility"
python "${REPO_ROOT}/extra_scripts/evaluate_suppa2_statistics.py" \
  --tealeaf-shards "${RUN_REPRO}/fold0/tealeaf_paired_path_uniform_eb_crossfit_total_shards" \
  --tealeaf-shards "${RUN_REPRO}/fold1/tealeaf_paired_path_uniform_eb_crossfit_total_shards" \
  --comparison "Tealeaf EC; SUPPA2 event definitions" \
    "${BASE}/fold0_merged/paired_path.tsv" \
    "${BASE}/fold1_merged/paired_path.tsv" \
  --output-dir "${BASE}/evaluation"
