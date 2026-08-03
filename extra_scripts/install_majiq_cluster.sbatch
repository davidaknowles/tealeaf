#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail
: "${MAJIQ_SOURCE:?set MAJIQ_SOURCE to the academic source checkout}"
MAJIQ_VENV=${MAJIQ_VENV:-${HOME}/venv/majiq_cluster}

module purge
module load Python/3.12.3-GCCcore-13.3.0 HTSlib/1.21 zstd/1.5.6-GCCcore-13.3.0 CMake/3.29.3-GCCcore-13.3.0
python -m venv "${MAJIQ_VENV}"
source "${MAJIQ_VENV}/bin/activate"
python -m pip install --upgrade pip
python -m pip install pybind11==2.13.6 scikit-build-core==0.11.1 ninja
htslib_root=$(dirname "$(dirname "$(command -v bgzip)")")
export HTSLIB_LIBRARY_DIR=${htslib_root}/lib
export HTSLIB_INCLUDE_DIR=${htslib_root}/include
python -m pip install --no-build-isolation "${MAJIQ_SOURCE}/moccasin"
python -m pip install --no-build-isolation "${MAJIQ_SOURCE}/majiq"
python -m pip install --no-build-isolation "${MAJIQ_SOURCE}/voila"
majiq-v3 --version
