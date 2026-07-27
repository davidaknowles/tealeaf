#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"

: "${ANALYSIS_TAG:?set ANALYSIS_TAG}"
: "${RANK:?set RANK}"
: "${INITIAL_FACTORS:?set INITIAL_FACTORS}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

GLM_RUN="${RUN_ROOT}/glm"
ALEVIN_DIR="${DATA_ROOT}/processed/merged_alevin"
GENE_LOSS_WEIGHT=${GENE_LOSS_WEIGHT:-0}

python -m tealeaf.sc.tealeaf_sc \
  --alevin_dir "${ALEVIN_DIR}" \
  --salmon_ref "${SALMON_REF}/spliceu.fa" \
  --ref_dir "${GLM_RUN}" \
  --cell_mode single_cell \
  --quant_method factorized \
  --regularization_target theta \
  --ec_design positional \
  --primer_pairs "${ALEVIN_DIR}/primer_pairs.tsv" \
  --min_half_umis 500 \
  --primer_sampling_model oligodt_tpm \
  --transcript_to_gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --gene_loss_weight "${GENE_LOSS_WEIGHT}" \
  --glm_device cuda \
  --glm_rank "${RANK}" \
  --glm_batch_cells 16384 \
  --glm_data_backend auto \
  --glm_no_write_chunks \
  --glm_no_minibatch \
  --glm_exact_inner_steps 32 \
  --glm_initial_factors "${INITIAL_FACTORS}" \
  --min_cell_umis 500 \
  --nucnorm_max_iter 8192 \
  --nucnorm_min_iter 50 \
  --nucnorm_patience 10 \
  --nucnorm_tol 1e-5 \
  --min_eq 5 \
  -o "${GLM_RUN}/${ANALYSIS_TAG}_theta_"
