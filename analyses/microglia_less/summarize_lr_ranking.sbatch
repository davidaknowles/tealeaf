#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=lr_rank_sensitivity
#SBATCH --output=.cache/lr_ranking_20261006/summary_%j.out
#SBATCH --error=.cache/lr_ranking_20261006/summary_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
pair_strength=$(python -c 'import json,sys; print(format(json.load(open(sys.argv[1]))["strength"], "g"))' "${RUN_ROOT}/differential/local_path_full_testing_eb_pairwise/selection/selection.json")
omnibus_strength=$(python -c 'import json,sys; print(format(json.load(open(sys.argv[1]))["strength"], "g"))' "${RUN_ROOT}/differential/local_path_full_testing_eb_omnibus_production/selection/selection.json")
python "${REPO_ROOT}/extra_scripts/summarize_tealeaf_lr_ranking.py" \
  --replication "${REPO_ROOT}/.cache/lr_ranking_20261006/pairwise_replication_precise.tsv" \
  --pairwise-tests "${RUN_ROOT}/differential/local_path_full_testing_eb_pairwise/a${pair_strength}/paired_path.tsv" \
  --omnibus-tests "${RUN_ROOT}/differential/local_path_full_testing_eb_omnibus_production/a${omnibus_strength}/paired_path.tsv" \
  --effects "${REPO_ROOT}/.cache/lr_ranking_20261006/effects" \
  --output-dir "${REPO_ROOT}/analyses/tilgner_long_read_replication/ranking_audit"
