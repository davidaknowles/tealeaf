#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=microglia_ec_block

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${SLURM_CPUS_PER_TASK}"

SHARD_COUNT=${SHARD_COUNT:-96}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
OUTPUT_TAG=${OUTPUT_TAG:-ec_block_glmm}
METHODS=${METHODS:-multinomial_full;multinomial_noise_full;dirichlet_multinomial_full}
NULL_REPLICATES=${NULL_REPLICATES:-20}
MIN_GENE_UMIS=${MIN_GENE_UMIS:-10}
MIN_GENE_SAMPLES=${MIN_GENE_SAMPLES:-0}
MIN_CELL_TYPES=${MIN_CELL_TYPES:-2}
MIN_CELLTYPE_MICE=${MIN_CELLTYPE_MICE:-3}
MAX_ITER=${MAX_ITER:-300}
RETRIES=${RETRIES:-2}
NULL_REPLICATE_RETRIES=${NULL_REPLICATE_RETRIES:-${RETRIES}}
VI_SAMPLES=${VI_SAMPLES:-16}
IFS=';' read -r -a method_args <<< "${METHODS}"

ARGS=(
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl"
  --features "${RUN_ROOT}/hierarchical/weighted_rank32_joint_glm_cols.txt"
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz"
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${SHARD_INDEX}"
  --methods "${method_args[@]}"
  --null-replicates "${NULL_REPLICATES}"
  --min-gene-umis "${MIN_GENE_UMIS}"
  --min-gene-samples "${MIN_GENE_SAMPLES}"
  --min-cell-types "${MIN_CELL_TYPES}"
  --min-celltype-mice "${MIN_CELLTYPE_MICE}"
  --max-iter "${MAX_ITER}"
  --retries "${RETRIES}"
  --null-replicate-retries "${NULL_REPLICATE_RETRIES}"
  --vi-samples "${VI_SAMPLES}"
  --shard-index "${SHARD_INDEX}"
  --shard-count "${SHARD_COUNT}"
)
if [[ -n "${EVENT_TABLE:-}" ]]; then
  ARGS+=(--event-table "${EVENT_TABLE}")
fi

python "${REPO_ROOT}/extra_scripts/run_ec_block_glmm.py" "${ARGS[@]}"
