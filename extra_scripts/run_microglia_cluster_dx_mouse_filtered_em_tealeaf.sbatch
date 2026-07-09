#!/bin/bash -l
#SBATCH -p bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH -t 24:00:00
#SBATCH -J tealeaf_microglia_em
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_min20_umi1e5_em_slurm-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_min20_umi1e5_em_slurm-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf

PREFIX=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_min20_umi1e5_em_
ALEVIN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/alevin_t2t_cache
BARCODES=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/barcodes_to_cluster_dx_mouse.filtered_min20_umi1e5.csv

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
  --thread 8
