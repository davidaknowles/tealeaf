#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 2-00:00:00
#SBATCH -J tealeaf_sc_glm
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

METHOD=${METHOD:?Set METHOD to admm_factorized, frank_wolfe, or factorized}
REGULARIZATION_TARGET=${REGULARIZATION_TARGET:-phi}
EC_DESIGN=${EC_DESIGN:-legacy}
ALEVIN=${ALEVIN:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin}
REF=${REF:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa}
REF_DIR=${REF_DIR:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref}
OUT_ROOT=${OUT_ROOT:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run}
GLM_MAX_ITER=${GLM_MAX_ITER:-2048}
FW_MAX_ATOMS=${FW_MAX_ATOMS:-256}
FW_POWER_ITER=${FW_POWER_ITER:-10}
FW_NONNEGATIVE_PENALTY=${FW_NONNEGATIVE_PENALTY:-1.0}

python -m tealeaf.sc.tealeaf_sc \
  --alevin_dir "${ALEVIN}" \
  --salmon_ref "${REF}" \
  --ref_dir "${REF_DIR}" \
  --cell_mode single_cell \
  --quant_method "${METHOD}" \
  --regularization_target "${REGULARIZATION_TARGET}" \
  --ec_design "${EC_DESIGN}" \
  --glm_device cuda \
  --glm_rank 64 \
  --fw_max_atoms "${FW_MAX_ATOMS}" \
  --fw_power_iter "${FW_POWER_ITER}" \
  --fw_nonnegative_penalty "${FW_NONNEGATIVE_PENALTY}" \
  --glm_batch_cells 4096 \
  --glm_no_write_chunks \
  --nucnorm_max_iter "${GLM_MAX_ITER}" \
  --nucnorm_min_iter 50 \
  --nucnorm_patience 10 \
  --nucnorm_tol 1e-6 \
  --min_eq 5 \
  -o "${OUT_ROOT}/single_cell_${EC_DESIGN}_${METHOD}_${REGULARIZATION_TARGET}_"
