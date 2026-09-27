#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=eval_primer_suppa
#SBATCH --output=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/primer-eval-%j.out
#SBATCH --error=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/primer-eval-%j.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

BASE="${REPO_ROOT}/.cache/suppa2_statistics"
for fold in 0 1; do
  python "${REPO_ROOT}/extra_scripts/filter_suppa2_primer_aware.py" \
    --input "${BASE}/primer_aware/all_tests.tsv.gz" \
    --fold "${fold}" \
    --output "${BASE}/primer_aware/filtered/fold${fold}.tsv.gz"
done

RUN_REPRO="${RUN_ROOT}/junction_benchmark/reproducibility"
python "${REPO_ROOT}/extra_scripts/evaluate_suppa2_statistics.py" \
  --tealeaf-shards "${RUN_REPRO}/fold0/tealeaf_paired_path_uniform_eb_crossfit_total_shards" \
  --tealeaf-shards "${RUN_REPRO}/fold1/tealeaf_paired_path_uniform_eb_crossfit_total_shards" \
  --comparison "SUPPA2 primer-aware; prior Wilcoxon" \
    "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2_primer_aware/split_data_fold0_tests.tsv.gz" \
    "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2_primer_aware/split_data_fold1_tests.tsv.gz" \
  --comparison "SUPPA2 primer-aware; alternative tests" \
    "${BASE}/primer_aware/filtered/fold0.tsv.gz" \
    "${BASE}/primer_aware/filtered/fold1.tsv.gz" \
  --output-dir "${BASE}/primer_aware/evaluation_alternatives"
