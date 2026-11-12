#!/bin/bash
# Within-subject label-swap null on GSE320507 (21 pairs of the seven most abundant cell types): LeafCutter (min_coverage 5)
# and Tealeaf (internal and terminal, settings of the pbmc_GSE320507_*_sh fits) on the same swaps.
set -euo pipefail
cd ~/projects/tealeaf
O=.cache/label_swap; PB=~/tealeaf_runs/pbmc/GSE320507/junction_benchmark
mkdir -p $O/logs
sbatch --parsable -p cpu --array=0-20 -J swap_lc -o $O/logs/lc.%A_%a.out --export="ALL,MIN_COVERAGE=5,REPO_ROOT=$PWD,BUNDLE=$PB/pseudobulk_junctions,CONTRASTS=$O/contrasts_swapped.json,COUNTS=$PB/leafcutter/clustering/benchmark_perind_numers.counts.gz,OUTPUT_DIR=$O/leafcutter,LEAFCUTTER_ROOT=$HOME/projects/leafcutter,PYTHON_MODULE=Python/3.12.3-GCCcore-13.3.0" extra_scripts/run_leafcutter_contrast.sbatch
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1'
COMMON="--subject-folds $PB/reproducibility/subject_folds.tsv --cohort full --level-pairs $O/level_pairs.tsv --swapped-pairs $O/swapped.tsv --min-share .01 --primer-offset --stored-opportunities --primer-offset-prior-sd 1 --shard-count 128"
for v in "internal|.cache/pbmc/GSE320507/collection/all|--structure internal" "terminal|.cache/pbmc/GSE320507/collection/junction|--structure terminal --group-same-junctions"; do
  IFS='|' read -r name reads flags <<< "$v"
  A=$(sbatch --parsable -p cpu --mem=16G -t 24:00:00 -a 0-127 -J swap_$name -o $O/logs/${name}_%A_%a.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/run_local_path_tests.py --local-reads $reads $flags $COMMON --output-dir $O/tealeaf_$name --shard-index \$SLURM_ARRAY_TASK_ID'")
  sbatch --parsable --dependency=afterok:$A -p cpu --mem=16G -t 1:00:00 -J swap_${name}_m -o $O/logs/${name}_merge_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort full --output-dir $O/tealeaf_$name --shard-count 128 --merge'"
  echo "$name $A"
done
