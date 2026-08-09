#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=16:00:00
#SBATCH --job-name=gse233208_ec_block
#SBATCH --output=/scratch/daknowles/gse233208_ec_block_%A_%a.out
#SBATCH --error=/scratch/daknowles/gse233208_ec_block_%A_%a.err

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

SHARD_COUNT=${SHARD_COUNT:-96}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
TEST_EFFECT=${TEST_EFFECT:-cell_type}
OUTPUT_TAG=${OUTPUT_TAG:-ec_block_${TEST_EFFECT}}
CANDIDATE_CACHE=${CANDIDATE_CACHE:?set CANDIDATE_CACHE}
CALIBRATION=${CALIBRATION:-bootstrap}
if [[ "${CALIBRATION}" == "lrt_bic" ]]; then
  DEFAULT_METHODS="laplace_multinomial;laplace_multinomial_noise"
else
  DEFAULT_METHODS="multinomial_full;multinomial_noise_full"
fi
IFS=';' read -r -a methods <<< "${METHODS:-${DEFAULT_METHODS}}"
ARGS=(
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl"
  --features "${RUN_ROOT}/hierarchical/full_rank32_joint_glm_cols.txt"
  --block-cache "${RUN_ROOT}/differential/gencode_v30_splice_blocks.json.gz"
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${SHARD_INDEX}"
  --methods "${methods[@]}" --test-effect "${TEST_EFFECT}"
  --calibration "${CALIBRATION}"
  --candidate-cache "${CANDIDATE_CACHE}"
  --null-replicates "${NULL_REPLICATES:-20}" --min-gene-umis "${MIN_GENE_UMIS:-10}"
  --min-celltype-mice "${MIN_CELLTYPE_SUBJECTS:-4}"
  --min-condition-mice "${MIN_CONDITION_SUBJECTS:-4}"
  --max-iter "${MAX_ITER:-300}" --retries "${RETRIES:-2}"
  --laplace-mode-steps "${LAPLACE_MODE_STEPS:-12}"
  --laplace-mode-gradient "${LAPLACE_MODE_GRADIENT:-implicit}"
  --laplace-cache-scope "${LAPLACE_CACHE_SCOPE:-gene_shape}"
  --latent-space "${LATENT_SPACE:-isoform}"
  --null-replicate-retries "${NULL_REPLICATE_RETRIES:-2}"
  --vi-samples "${VI_SAMPLES:-16}"
  --shard-index "${SHARD_INDEX}" --shard-count "${SHARD_COUNT}"
)
if [[ -n "${MAX_CANDIDATES:-}" ]]; then
  ARGS+=(--max-candidates "${MAX_CANDIDATES}")
fi
if [[ -n "${SUBJECT_FOLDS:-}" ]]; then
  ARGS+=(--subject-folds "${SUBJECT_FOLDS}" --subject-fold "${SUBJECT_FOLD:?}")
fi
python "${REPO_ROOT}/extra_scripts/run_ec_block_glmm.py" "${ARGS[@]}"
