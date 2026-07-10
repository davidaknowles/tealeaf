#!/bin/bash -l
#SBATCH -p bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH -t 2-00:00:00
#SBATCH -J tealeaf_weighted_glm
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/weighted_all_glm_slurm-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/weighted_all_glm_slurm-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf

ALEVIN=${ALEVIN:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin}
BARCODES=${BARCODES:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/barcodes_to_cluster_dx_mouse.filtered_min20_umi1e5.csv}
OUT_ROOT=${OUT_ROOT:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run}
METHODS=${METHODS:-"em nnls nnls_nucnorm"}

for METHOD in ${METHODS}; do
  PREFIX="${OUT_ROOT}/weighted_parsimony_min20_umi1e5_${METHOD}_"
  python /gpfs/commons/home/daknowles/projects/tealeaf/tealeaf/sc/tealeaf_sc.py \
    --alevin_dir "${ALEVIN}" \
    --salmon_ref /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa \
    --ref_dir /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref \
    --ref_prefix vM32_ \
    --pseudobulk_samples "${BARCODES}" \
    -o "${PREFIX}" \
    --min_eq 5 \
    --normalization_scale junction \
    --introncutoff 80 \
    --minclucounts 100 \
    --thread 8 \
    --quant_method "${METHOD}"
done
