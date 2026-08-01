#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-47
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=microglia_ec_glmm

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${SLURM_CPUS_PER_TASK}"

SHARD_COUNT=${SHARD_COUNT:-48}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
CONTRAST=${CONTRAST:-condition}
OUTPUT_TAG=${OUTPUT_TAG:-ec_glmm_${CONTRAST}}
METHODS=${METHODS:-laplace_multinomial;tilted_elbo;renyi_multinomial;laplace_dirichlet_multinomial;elbo_dirichlet_multinomial;renyi_dirichlet_multinomial}
MAX_ITER=${MAX_ITER:-200}
VI_SAMPLES=${VI_SAMPLES:-128}
EVALUATION_SAMPLES=${EVALUATION_SAMPLES:-2048}
MAX_GENES=${MAX_GENES:-}
IFS=';' read -r -a method_args <<< "${METHODS}"
max_genes_args=()
if [[ -n "${MAX_GENES}" ]]; then
  max_genes_args+=(--max-genes "${MAX_GENES}")
fi

python "${REPO_ROOT}/extra_scripts/run_ec_glmm.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${RUN_ROOT}/primer_pairs.tsv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --barcode-groups "${RUN_ROOT}/barcodes_to_cluster_dx_mouse.filtered_min20_umi1e5.csv" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${SHARD_INDEX}" \
  --contrast "${CONTRAST}" \
  --methods "${method_args[@]}" \
  --max-iter "${MAX_ITER}" \
  --vi-samples "${VI_SAMPLES}" \
  --evaluation-samples "${EVALUATION_SAMPLES}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}" \
  "${max_genes_args[@]}"
