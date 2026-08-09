#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --job-name=gse233208_path_wald
#SBATCH --output=/scratch/daknowles/gse233208_path_wald_%A_%a.out
#SBATCH --error=/scratch/daknowles/gse233208_path_wald_%A_%a.err

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
: "${CANDIDATE_CACHE:?set CANDIDATE_CACHE}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

SHARD_COUNT=${SHARD_COUNT:-96}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
OUTPUT_TAG=${OUTPUT_TAG:-path_wald}

python "${REPO_ROOT}/extra_scripts/run_path_wald_ablation.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE}" \
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${SHARD_INDEX}" \
  --covariance "${COVARIANCE:-conditional}" \
  --fit-max-iter "${FIT_MAX_ITER:-100}" \
  --label-permutations "${LABEL_PERMUTATIONS:-0}" \
  --wild-null-replicates "${WILD_NULL_REPLICATES:-0}" \
  --seed "${SEED:-0}" \
  --inference "${INFERENCE:-wald}" \
  --effect-prior-scales "${EFFECT_PRIOR_SCALES:-0.1;0.25;0.5;1.0}" \
  --gaussian-solver "${GAUSSIAN_SOLVER:-woodbury}" \
  --cluster-adjustment "${CLUSTER_ADJUSTMENT:-cr1}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}"
