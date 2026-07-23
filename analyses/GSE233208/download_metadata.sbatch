#!/bin/bash
#SBATCH --job-name=gse233208_metadata
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00

set -euo pipefail
source "$(dirname "$0")/config.env"
destination="${DATA_ROOT}/metadata"
mkdir -p "${destination}"
base=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/suppl
for filename in \
  GSE233208_AD-DS_Cases.csv.gz \
  GSE233208_Human_snRNA-Seq_ADDS_integrated.rds.gz
do
  curl --fail --location --retry 8 --retry-all-errors \
    --continue-at - \
    --output "${destination}/${filename}.part" \
    "${base}/${filename}"
  mv "${destination}/${filename}.part" "${destination}/${filename}"
done
