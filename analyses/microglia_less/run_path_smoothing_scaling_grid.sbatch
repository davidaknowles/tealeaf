#!/bin/bash -l
#SBATCH -p dev
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=path_smoothing_scaling

set -euo pipefail

source "${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export REPO_ROOT
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export ALPHA_GRID="${ALPHA_GRID:-0.5;1;2;4;8;16;32;64;128}"
export MAX_CANDIDATES="${MAX_CANDIDATES:-1000}"
export SHARD_COUNT="${SHARD_COUNT:-96}"

repro_root="${RUN_ROOT}/junction_benchmark/reproducibility"
output_tag="${OUTPUT_TAG:-path_smoothing_eb3}"
workers="${SLURM_CPUS_PER_TASK:-1}"

for scaling in per_path total; do
  export PATH_PSEUDOCOUNT_SCALING="${scaling}"
  for fold in 0 1; do
    fold_root="${repro_root}/fold${fold}"
    export CANDIDATE_CACHE="${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold${fold}_candidates.pkl"
    export OUTPUT_ROOT="${fold_root}/${output_tag}_${scaling}_evidence"
    mkdir -p "${OUTPUT_ROOT}"
    seq 0 $((SHARD_COUNT - 1)) | xargs -P "${workers}" -I '{}' bash -c 'SLURM_ARRAY_TASK_ID="$1" bash "$REPO_ROOT/analyses/microglia_less/run_path_smoothing.sbatch"' _ '{}'
    python "${REPO_ROOT}/extra_scripts/merge_path_smoothing.py" --inputs "${OUTPUT_ROOT}"/evidence_*.tsv --output "${fold_root}/${output_tag}_${scaling}.json"
  done
done
