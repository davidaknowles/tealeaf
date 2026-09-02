#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=full_path_eb_select

set -euo pipefail

source "${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

IFS=';' read -r -a strengths <<< "${STRENGTHS:-2;4;8;16;32;64;128;256}"
args=()
for strength in "${strengths[@]}"; do
  args+=(
    --result
    "${strength}"
    "${RESULT_ROOT:?}/a${strength}"
  )
done

python "${REPO_ROOT}/extra_scripts/select_path_smoothing_from_pvalues.py" \
  "${args[@]}" \
  --output-dir "${RESULT_ROOT}/selection"

python - "${RESULT_ROOT}" <<'PY'
import json
import pathlib
import shutil
import sys

root = pathlib.Path(sys.argv[1])
selection = json.loads((root / "selection/selection.json").read_text())
strength = selection["strength"]
label = str(int(strength)) if float(strength).is_integer() else str(strength)
shutil.copy2(root / f"a{label}/summary.tsv", root / "selection/selected_summary.tsv")
PY
