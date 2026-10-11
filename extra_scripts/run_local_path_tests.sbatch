#!/bin/bash
#SBATCH --job-name=local_path_tests
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --array=0-95
#SBATCH --output=.cache/binary_grid/logs/local_tests_%A_%a.log
# Task t runs shard t%32 of cohort fold0/fold1/full. MODEL names the output;
# EXTRA_ARGS selects the variant. Folds test the matched split universes; full
# data also tests every pair of long-read-mapped cell types.
set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/cpa_blackwell/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
folds=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/junction_benchmark/reproducibility/subject_folds.tsv
cohorts=(0 1 full)
cohort=${cohorts[$((SLURM_ARRAY_TASK_ID / 32))]}
shard=$((SLURM_ARRAY_TASK_ID % 32))
lr=()
if [ "$cohort" = full ]; then lr=(--lr-cell-types "ASC;EX_cerebellum_granule;EX_cortical;EX_hippocampus_granule;EX_hippocampus_pyramidal;EX_thalamus;INH_medium_spiny;ODC"); fi
python -W ignore extra_scripts/run_local_path_tests.py --local-reads "${LOCAL_READS:-.cache/local_path_reads_annotated}" --subject-folds "$folds" --cohort "$cohort" --output-dir ".cache/binary_grid/${MODEL:?}_${cohort}/shard_${shard}" --shard-index "$shard" --shard-count 32 --universe analyses/split_coverage_direction/event_coverage_pvalues.tsv.gz "${lr[@]}" ${EXTRA_ARGS:-}
