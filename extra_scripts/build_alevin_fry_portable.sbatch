#!/bin/bash -l
#SBATCH --job-name=alevin_fry_portable_build
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/build_logs/alevin_fry_portable_build-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/build_logs/alevin_fry_portable_build-%j.err

set -euo pipefail

module load CMake/3.26.3-GCCcore-12.3.0 GCC/12.3.0

cd /gpfs/commons/home/daknowles/projects/alevin-fry

cargo build \
  --release \
  --target-dir target-portable \
  --config "build.rustflags=[]" \
  --config .cargo/config-portable.toml
