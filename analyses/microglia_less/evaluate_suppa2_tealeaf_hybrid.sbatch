#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=eval_suppa_tealeaf
#SBATCH --output=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/hybrid-eval-%j.out
#SBATCH --error=/scratch/daknowles/tealeaf-suppa2-statistics-20260927/hybrid-eval-%j.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

RUN_REPRO="${RUN_ROOT}/junction_benchmark/reproducibility"
HYBRID="${REPO_ROOT}/.cache/suppa2_statistics/tealeaf_hybrid_gene_wide"
python "${REPO_ROOT}/extra_scripts/evaluate_suppa2_statistics.py" \
  --tealeaf-tests "${RUN_REPRO}/fold0/tealeaf_paired_path_total_a32_production/paired_path.tsv" \
  --tealeaf-tests "${RUN_REPRO}/fold1/tealeaf_paired_path_total_a32_production/paired_path.tsv" \
  --comparison "SUPPA2 definitions; Tealeaf EC model" \
    "${HYBRID}/fold0_merged/paired_path.tsv" \
    "${HYBRID}/fold1_merged/paired_path.tsv" \
  --comparison "Tealeaf direct local path" \
    "${RUN_REPRO}/fold0/tealeaf_paired_path_total_a32_production/paired_path.tsv" \
    "${RUN_REPRO}/fold1/tealeaf_paired_path_total_a32_production/paired_path.tsv" \
  --comparison "SUPPA2 native; all event types" \
    "${REPO_ROOT}/.cache/suppa2_statistics/native_all_events/fold0.tsv.gz" \
    "${REPO_ROOT}/.cache/suppa2_statistics/native_all_events/fold1.tsv.gz" \
  --output-dir "${HYBRID}/evaluation_expanded"
