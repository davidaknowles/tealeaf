#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=microglia_ec_block
#SBATCH --output=/scratch/daknowles/tealeaf_ec_block_%A_%a.out
#SBATCH --error=/scratch/daknowles/tealeaf_ec_block_%A_%a.err

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

unset LOADEDMODULES
for variable in $(compgen -A variable _ModuleTable); do
  unset "${variable}"
done
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${SLURM_CPUS_PER_TASK}"

SHARD_COUNT=${SHARD_COUNT:-96}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
OUTPUT_TAG=${OUTPUT_TAG:-ec_block_glmm}
CALIBRATION=${CALIBRATION:-bootstrap}
if [[ "${CALIBRATION}" == "lrt_bic" ]]; then
  DEFAULT_METHODS="laplace_multinomial;laplace_multinomial_noise"
else
  DEFAULT_METHODS="multinomial_full;multinomial_noise_full;dirichlet_multinomial_full"
fi
METHODS=${METHODS:-${DEFAULT_METHODS}}
NULL_REPLICATES=${NULL_REPLICATES:-20}
TEST_EFFECT=${TEST_EFFECT:-cell_type}
MIN_GENE_UMIS=${MIN_GENE_UMIS:-10}
MIN_GENE_SAMPLES=${MIN_GENE_SAMPLES:-0}
MIN_CELL_TYPES=${MIN_CELL_TYPES:-2}
MIN_CELLTYPE_MICE=${MIN_CELLTYPE_MICE:-3}
MIN_CONDITIONS=${MIN_CONDITIONS:-2}
MIN_CONDITION_MICE=${MIN_CONDITION_MICE:-3}
MAX_ITER=${MAX_ITER:-300}
LAPLACE_MODE_STEPS=${LAPLACE_MODE_STEPS:-12}
RETRIES=${RETRIES:-2}
NULL_REPLICATE_RETRIES=${NULL_REPLICATE_RETRIES:-${RETRIES}}
VI_SAMPLES=${VI_SAMPLES:-16}
CANDIDATE_CACHE=${CANDIDATE_CACHE:-}
PAIRWISE_NULL_SEED=${PAIRWISE_NULL_SEED:-}
MAX_CANDIDATES=${MAX_CANDIDATES:-}
JOINT_GENE_TEST=${JOINT_GENE_TEST:-0}
IFS=';' read -r -a method_args <<< "${METHODS}"

ARGS=(
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl"
  --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt"
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz"
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${SHARD_INDEX}"
  --methods "${method_args[@]}"
  --calibration "${CALIBRATION}"
  --null-replicates "${NULL_REPLICATES}"
  --test-effect "${TEST_EFFECT}"
  --min-gene-umis "${MIN_GENE_UMIS}"
  --min-gene-samples "${MIN_GENE_SAMPLES}"
  --min-cell-types "${MIN_CELL_TYPES}"
  --min-celltype-mice "${MIN_CELLTYPE_MICE}"
  --min-conditions "${MIN_CONDITIONS}"
  --min-condition-mice "${MIN_CONDITION_MICE}"
  --max-iter "${MAX_ITER}"
  --laplace-mode-steps "${LAPLACE_MODE_STEPS}"
  --laplace-mode-gradient "${LAPLACE_MODE_GRADIENT:-implicit}"
  --laplace-cache-scope "${LAPLACE_CACHE_SCOPE:-gene_shape}"
  --laplace-optimizer "${LAPLACE_OPTIMIZER:-lbfgs}"
  --adam-learning-rate "${ADAM_LEARNING_RATE:-0.03}"
  --adam-steps "${ADAM_STEPS:-100}"
  --null-warm-start "${NULL_WARM_START:-default}"
  --latent-space "${LATENT_SPACE:-isoform}"
  --retries "${RETRIES}"
  --null-replicate-retries "${NULL_REPLICATE_RETRIES}"
  --vi-samples "${VI_SAMPLES}"
  --shard-index "${SHARD_INDEX}"
  --shard-count "${SHARD_COUNT}"
)
if [[ -n "${EVENT_TABLE:-}" ]]; then
  ARGS+=(--event-table "${EVENT_TABLE}")
fi
if [[ -n "${CANDIDATE_CACHE}" ]]; then
  ARGS+=(--candidate-cache "${CANDIDATE_CACHE}")
fi
if [[ -n "${SUBJECT_FOLDS:-}" ]]; then
  ARGS+=(--subject-folds "${SUBJECT_FOLDS}" --subject-fold "${SUBJECT_FOLD:?}")
fi
if [[ -n "${PAIRWISE_NULL_SEED}" ]]; then
  ARGS+=(--pairwise-null-seed "${PAIRWISE_NULL_SEED}")
fi
if [[ -n "${MAX_CANDIDATES}" ]]; then
  ARGS+=(--max-candidates "${MAX_CANDIDATES}")
fi
if [[ "${JOINT_GENE_TEST}" == "1" ]]; then
  ARGS+=(--joint-gene-test)
fi

python "${REPO_ROOT}/extra_scripts/run_ec_block_glmm.py" "${ARGS[@]}"
