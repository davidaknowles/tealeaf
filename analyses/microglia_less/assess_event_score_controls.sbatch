#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --job-name=assess_event_score
#SBATCH --output=.cache/event_information_20261007/assess_event_score_%j.out

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
INPUT_ROOT=${EVENT_SCORE_INPUT_ROOT:-${REPO_ROOT}/.cache/event_information_20261007}
INFORMATION_METRIC=${EVENT_SCORE_INFORMATION_METRIC:-absolute}
OUTPUT=${EVENT_SCORE_OUTPUT:-${REPO_ROOT}/analyses/event_mixed_score}
python "${REPO_ROOT}/extra_scripts/assess_event_input_controls.py" \
  --cache "${INPUT_ROOT}" --source original_binary --output-dir "${OUTPUT}" \
  --matrix-dir "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" \
  --gtf "${GENCODE_GTF}" \
  --event-catalog "${REPO_ROOT}/.cache/suppa2_statistics/all_event_catalog/event_catalog.tsv.gz" \
  --inference mixed-score --shard-count 64 --information-metric "${INFORMATION_METRIC}"
