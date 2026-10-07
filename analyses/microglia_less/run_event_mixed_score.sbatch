#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-191%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --job-name=event_mixed_score
#SBATCH --output=.cache/event_information_20261007/mixed_%A_%a.out
#SBATCH --error=.cache/event_information_20261007/mixed_%A_%a.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SHARD_COUNT=64
COHORT_INDEX=$((SLURM_ARRAY_TASK_ID / SHARD_COUNT))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARD_COUNT))
COHORTS=(full fold0 fold1)
COHORT=${COHORTS[$COHORT_INDEX]}
if [[ "${COHORT}" == "full" ]]; then
  CANDIDATES="${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl"
else
  CANDIDATES="${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_${COHORT}_candidates.pkl"
fi
CONTROL="${REPO_ROOT}/.cache/event_information_20261007/original_binary_paired"
python "${REPO_ROOT}/extra_scripts/run_suppa2_tealeaf_hybrid.py" \
  --data-cache "${CONTROL}/prepared.pkl" --features "${CONTROL}/features.txt" \
  --candidate-cache "${CANDIDATES}" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --output-dir "${REPO_ROOT}/.cache/event_information_20261007/mixed_score/original_binary/${COHORT}/shard_${SHARD_INDEX}" \
  --scan-all-events --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}" \
  --inference mixed-score --path-pseudocount 0 --max-iter 2000 --null-multistart \
  --profile-event-mass --report-pseudocount 1 --null-replicates 32
