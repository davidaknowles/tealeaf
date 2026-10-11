#!/bin/bash
#SBATCH --job-name=assess_local
#SBATCH --partition=cpu
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --output=.cache/binary_grid/logs/assess_local_%j.log
# Merge the three cohorts of a run_local_path_tests.sbatch model and run the
# split and long-read assessments into analyses/$MODEL.
set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/cpa_blackwell/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1
R=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/differential
for c in 0 1 full; do
  python extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort "$c" --output-dir ".cache/binary_grid/${MODEL:?}_$c" --shard-count 32 --merge
done
python -W ignore extra_scripts/assess_paired_inference_audit.py --cache .cache/binary_grid --model "$MODEL" --output-dir "analyses/$MODEL" --matrix-dir .cache/tilgner_scisorseq/LR_annotTranscripts_byCT --gtf /gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor/gencode.vM32.basic.annotation.gtf --block-cache "$R/gencode_vM32_splice_blocks.json.gz" --skip-gene-umi-check
