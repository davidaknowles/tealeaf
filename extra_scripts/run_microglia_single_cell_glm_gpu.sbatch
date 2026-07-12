#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 2-00:00:00
#SBATCH -J tealeaf_sc_glm

set -euo pipefail

module load PyTorch/2.1.2-foss-2023b SciPy-bundle/2023.11-gfbf-2023b
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

METHOD=${METHOD:?Set METHOD to admm_factorized, frank_wolfe, or factorized}
REGULARIZATION_TARGET=${REGULARIZATION_TARGET:-phi}
ALEVIN=${ALEVIN:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin}
REF=${REF:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa}
REF_DIR=${REF_DIR:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref}
OUT_ROOT=${OUT_ROOT:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run}

python -m tealeaf.sc.tealeaf_sc \
  --alevin_dir "${ALEVIN}" \
  --salmon_ref "${REF}" \
  --ref_dir "${REF_DIR}" \
  --cell_mode single_cell \
  --quant_method "${METHOD}" \
  --regularization_target "${REGULARIZATION_TARGET}" \
  --glm_device cuda \
  --glm_rank 64 \
  --glm_batch_cells 4096 \
  --nucnorm_max_iter 64 \
  --nucnorm_tol 1e-4 \
  --min_eq 5 \
  -o "${OUT_ROOT}/single_cell_weighted_${METHOD}_${REGULARIZATION_TARGET}_"
