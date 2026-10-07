#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-15%16
#SBATCH --job-name=ec_count_null
#SBATCH --output=.cache/paired_reporting_20261006/count_null_%A_%a.out
#SBATCH --error=.cache/paired_reporting_20261006/count_null_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ -n "${RESIDUAL_CONCENTRATION:-}" ]]; then
  ARGS+=(--residual-concentration "${RESIDUAL_CONCENTRATION}")
fi
if [[ "${JOINT_DM_ONLY:-0}" == 1 ]]; then
  ARGS+=(--joint-dm-only)
fi
if [[ "${COX_REID_ONLY:-0}" == 1 ]]; then
  ARGS+=(--cox-reid-only)
fi
if [[ -n "${KNOWN_CONCENTRATION:-}" ]]; then
  ARGS+=(--known-concentration "${KNOWN_CONCENTRATION}")
fi
if [[ "${FREE_ISOFORMS:-0}" == 1 ]]; then
  ARGS+=(--free-isoforms)
fi
if [[ "${INCLUDE_SCORE:-0}" == 1 ]]; then
  ARGS+=(--include-score)
fi
python "${REPO_ROOT}/extra_scripts/audit_ec_count_null.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE:-${RUN_ROOT}/differential/ec_block_reproducibility_fold0_candidates.pkl}" \
  --reference "${REFERENCE:-${REPO_ROOT}/analyses/path_reporting_omnibus/omnibus_0_tests.tsv.gz}" \
  --output-dir "${OUTPUT_ROOT:-${REPO_ROOT}/.cache/paired_reporting_20261006/ec_count_null}/shard_${SLURM_ARRAY_TASK_ID}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count 16 \
  --blocks "${NULL_BLOCKS:-128}" --draws "${NULL_DRAWS:-16}" --null-replicates 32 \
  --subject-scale "${SUBJECT_SCALE:-.5}" --mode "${MODE:-omnibus}" --prior-center "${PRIOR_CENTER:-uniform}" "${ARGS[@]}"
