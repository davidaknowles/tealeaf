#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=1-00:00:00

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load R/4.3.2
Rscript "${REPO_ROOT}/extra_scripts/export_seurat_gene_pca.R" \
  "${DATA_ROOT}/metadata/GSE233208_Human_snRNA-Seq_ADDS_integrated.rds.gz" \
  "${DATA_ROOT}/metadata/GSE233208_Human_snRNA-Seq_ADDS_gene_pca.tsv.gz" \
  2000 \
  30
