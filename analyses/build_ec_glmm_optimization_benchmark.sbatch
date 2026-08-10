#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --job-name=ec_benchmark_build
#SBATCH --output=/gpfs/commons/home/daknowles/tealeaf_runs/ec_benchmark_build_%j.out

set -euo pipefail

: "${DATA_CACHE:?set DATA_CACHE}"
: "${CANDIDATE_CACHE:?set CANDIDATE_CACHE}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"

cd /gpfs/commons/home/daknowles/projects/tealeaf
module load Python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=.
/gpfs/commons/home/daknowles/venv/jax/bin/python analyses/build_ec_glmm_optimization_benchmark.py --data-cache "${DATA_CACHE}" --candidate-cache "${CANDIDATE_CACHE}" --output-dir "${OUTPUT_DIR}"
