#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=ec_analytic_test
#SBATCH --output=/gpfs/commons/home/daknowles/tealeaf_runs/ec_analytic_test_%j.out

set -euo pipefail

cd /gpfs/commons/home/daknowles/projects/tealeaf
module load Python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=.
export JAX_PLATFORMS=cpu
/gpfs/commons/home/daknowles/venv/jax/bin/python -m pytest -q ${PYTEST_TARGETS:-tests/test_ec_glmm.py tests/test_ec_block_glmm.py}
