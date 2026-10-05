#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=rmats_split_publish
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/publish_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/publish_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/publish_rmats_split_data.py" --repo-root "${REPO_ROOT}" --expected-document-sha "${EXPECTED_DOCUMENT_SHA}"
cd "${REPO_ROOT}/docs"
/gpfs/commons/home/daknowles/knowles_lab/software/texlive/2026/bin/x86_64-linux/latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../.cache/latex_coverage_audit differential.tex
cp ../.cache/latex_coverage_audit/differential.pdf differential.pdf
