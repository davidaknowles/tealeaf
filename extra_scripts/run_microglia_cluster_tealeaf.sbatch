#!/bin/bash -l
#SBATCH -p bigmem
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH -t 08:00:00
#SBATCH -J tealeaf_microglia
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_slurm-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_slurm-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf

PREFIX=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_
ALEVIN=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu/out_permit_known/quant_t2t_dedup/alevin

python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/aggregate_pseudobulk_sparse.py \
  --matrix "${ALEVIN}/quants_mat.npz" \
  --rows "${ALEVIN}/quants_mat_rows.txt" \
  --cols "${ALEVIN}/quants_mat_cols.txt" \
  --barcode-groups /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/barcodes_to_cluster_name.csv \
  -o "${PREFIX}"

python /gpfs/commons/home/daknowles/projects/tealeaf/tealeaf/sc/tealeaf_sc.py \
  --preprocessed \
  --alevin_dir "${ALEVIN}" \
  --salmon_ref /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa \
  --ref_dir /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref \
  --ref_prefix vM32_ \
  -o "${PREFIX}" \
  --introncutoff 80 \
  --minclucounts 100 \
  --thread 4
