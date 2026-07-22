#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 2-00:00:00
#SBATCH -J tealeaf_paired_glm
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

DESIGN=${DESIGN:?Set DESIGN to binary or weighted}
RANK=${RANK:-64}
RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run
ALEVIN=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin

python -m tealeaf.sc.tealeaf_sc \
  --alevin_dir "${ALEVIN}" \
  --salmon_ref /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa \
  --ref_dir /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref \
  --cell_mode single_cell \
  --primer_pairs "${RUN}/primer_pairs.tsv" \
  --min_half_umis 500 \
  --quant_method factorized \
  --regularization_target theta \
  --ec_design "${DESIGN}" \
  --glm_device cuda \
  --glm_rank "${RANK}" \
  --glm_batch_cells 16384 \
  --glm_data_backend auto \
  --glm_polish_max_iter 32 \
  --glm_no_write_chunks \
  --nucnorm_max_iter 2048 \
  --nucnorm_min_iter 50 \
  --nucnorm_patience 10 \
  --nucnorm_tol 1e-5 \
  --min_eq 5 \
  -o "${RUN}/paired_primer_${DESIGN}_factorized_theta_"
