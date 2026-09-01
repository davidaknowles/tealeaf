#!/bin/bash -l
#SBATCH -p dev
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=path_eb_scaling

set -euo pipefail

source "${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export REPO_ROOT
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NULL_REPLICATES="${NULL_REPLICATES:-32}"
export PATH_PRIOR_CENTER=baseline
export RETAIN_UNCERTAINTY=0
export SHARD_COUNT="${SHARD_COUNT:-96}"

benchmark_root="${RUN_ROOT}/junction_benchmark"
repro_root="${benchmark_root}/reproducibility"
output_tag="${OUTPUT_TAG:-tealeaf_paired_path_eb3}"
workers="${SLURM_CPUS_PER_TASK:-1}"

for scaling in per_path total; do
  export PATH_PSEUDOCOUNT_SCALING="${scaling}"
  for fold in 0 1; do
    fold_root="${repro_root}/fold${fold}"
    export CANDIDATE_CACHE="${RUN_ROOT}/differential/ec_block_reproducibility_pairwise_fold${fold}_candidates.pkl"
    export SMOOTHING_MAP="${fold_root}/path_smoothing_eb3_${scaling}.json"
    export OUTPUT_ROOT="${fold_root}/${output_tag}_${scaling}_shards"
    mkdir -p "${OUTPUT_ROOT}"
    seq 0 $((SHARD_COUNT - 1)) | xargs -P "${workers}" -I '{}' bash -c 'SLURM_ARRAY_TASK_ID="$1" bash "$REPO_ROOT/analyses/microglia_less/run_paired_path_test.sbatch"' _ '{}'
    python "${REPO_ROOT}/extra_scripts/merge_paired_path_test.py" --shards "${OUTPUT_ROOT}"/shard_* --output-dir "${fold_root}/${output_tag}_${scaling}" --moderate-variances --calibration empirical
  done
  python "${REPO_ROOT}/extra_scripts/score_ds_reproducibility.py" --fold-dir "${repro_root}/fold0" --fold-dir "${repro_root}/fold1" --junction-table "${benchmark_root}/pseudobulk_junctions_junctions.tsv.gz" --leafcutter-counts "${benchmark_root}/leafcutter/clustering/benchmark_perind_numers.counts.gz" --leafcutter-map "${repro_root}/leafcutter_cluster_gene.tsv.gz" --tealeaf-subdir "${output_tag}_${scaling}" --comparison-input-subdir comparison_majiq_min3_cov3 --tealeaf-label "Tealeaf EB ${scaling}" --tealeaf-format paired_path --min-paired-subjects 4 --match-pairs --output-dir "${repro_root}/${output_tag}_${scaling}_comparison"
done
