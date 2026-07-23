#!/bin/bash
#SBATCH --job-name=gse233208_export_meta
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load R/4.3.2
Rscript "${REPO_ROOT}/extra_scripts/export_seurat_metadata.R" \
  "${DATA_ROOT}/metadata/GSE233208_Human_snRNA-Seq_ADDS_integrated.rds.gz" \
  "${DATA_ROOT}/metadata/GSE233208_Human_snRNA-Seq_ADDS_metadata.tsv.gz"
