#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=hybrid_split_eval
#SBATCH --output=.cache/suppa2_statistics/hybrid-split-eval-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-split-eval-%j.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
BASE=${DEBUG_ROOT:-.cache/suppa2_statistics/hybrid_lr_analytic}/profile32_split
for FOLD in 0 1; do
  python extra_scripts/merge_paired_path_test.py --shards "${BASE}/fold${FOLD}"/shard_* --output-dir "${BASE}/fold${FOLD}_merged" --calibration empirical
done
RUN_REPRO="${RUN_ROOT}/junction_benchmark/reproducibility"
python extra_scripts/evaluate_suppa2_statistics.py \
  --tealeaf-shards "${RUN_REPRO}/fold0/tealeaf_paired_path_uniform_eb_crossfit_total_shards" \
  --tealeaf-shards "${RUN_REPRO}/fold1/tealeaf_paired_path_uniform_eb_crossfit_total_shards" \
  --comparison "Tealeaf EC; SUPPA2 event definitions, profiled event mass" \
    "${BASE}/fold0_merged/paired_path.tsv" "${BASE}/fold1_merged/paired_path.tsv" \
  --output-dir "${BASE}/evaluation"
