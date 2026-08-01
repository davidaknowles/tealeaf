#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-47
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --job-name=microglia_cluster_models

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${SLURM_CPUS_PER_TASK}"

SOURCE_TAG=${SOURCE_TAG:-weighted_hierarchical_subisoform_min25}
OUTPUT_TAG=${OUTPUT_TAG:-clustered_compositional_models}
SHARD_COUNT=${SHARD_COUNT:-48}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
PERMUTATIONS=${PERMUTATIONS:-50}
POWER_REPLICATES=${POWER_REPLICATES:-50}
MAX_ITER=${MAX_ITER:-200}

python "${REPO_ROOT}/extra_scripts/run_clustered_compositional_models.py" \
  --estimates "${RUN_ROOT}/differential/${SOURCE_TAG}/subisoform_estimates.jsonl.gz" \
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz" \
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${SHARD_INDEX}" \
  --min-condition-replicates 3 \
  --permutations "${PERMUTATIONS}" \
  --power-replicates "${POWER_REPLICATES}" \
  --effect-sizes 0.25 0.5 \
  --max-iter "${MAX_ITER}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}"
