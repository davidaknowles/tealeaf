#!/bin/bash
# Label-swap null on both halves of GSE320507 over all 171 cell-type pairs for one swap draw (make_contrasts_spread.py
# --all-pairs --seed S --output-dir O): LeafCutter (min_coverage 5), scQuint, rMATS (paired PAIRADISE JCEC) and Tealeaf
# (internal and terminal, settings of the pbmc_GSE320507_*_sh half fits), all with the same swaps.
# Usage: submit_all.sh <draw dir, e.g. .cache/label_swap_all/seed1>
set -euo pipefail
cd ~/projects/tealeaf
O=$1; PB=~/tealeaf_runs/pbmc/GSE320507/junction_benchmark; S=128; D=$(basename $O)
mkdir -p $O/logs
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1'
for k in 0 1; do
  C=$O/fold${k}_junction_contrasts.json
  sbatch --parsable -p cpu --array=0-170%60 -J sw_${D}_lc$k -o $O/logs/lc$k.%A_%a.out --export="ALL,MIN_COVERAGE=5,REPO_ROOT=$PWD,BUNDLE=$PB/pseudobulk_junctions,CONTRASTS=$C,COUNTS=$PB/leafcutter/clustering/benchmark_perind_numers.counts.gz,OUTPUT_DIR=$O/fold$k/leafcutter,LEAFCUTTER_ROOT=$HOME/projects/leafcutter,PYTHON_MODULE=Python/3.12.3-GCCcore-13.3.0" extra_scripts/run_leafcutter_contrast.sbatch
  sbatch --parsable -p cpu --cpus-per-task=4 --mem=32G -t 12:00:00 --array=0-170%60 -J sw_${D}_scq$k -o $O/logs/scq$k.%A_%a.out --export="ALL,REPO_ROOT=$PWD,BUNDLE=$PB/pseudobulk_junctions,CONTRASTS=$C,OUTPUT_DIR=$O/fold$k/scquint,SCQUINT_ROOT=$HOME/projects/scquint,PYTHON_MODULE=Python/3.12.3-GCCcore-13.3.0" extra_scripts/run_scquint_contrast.sbatch
  sbatch --parsable --array=0-170%80 -J sw_${D}_rmats$k -o $O/logs/rmats$k.%A_%a.out --export="ALL,CONTRASTS=$PWD/$O/fold${k}_rmats_contrasts.json,OUTDIR=$PWD/$O/fold$k/rmats" analyses/label_swap_null/rmats_swap.sbatch
  COMMON="--subject-folds .cache/pbmc/GSE320507/subject_folds.tsv --cohort $k --swapped-pairs $O/swapped.tsv --min-share .01 --primer-offset --stored-opportunities --primer-offset-prior-sd 1 --shard-count $S"
  for v in "internal|.cache/pbmc/GSE320507/collection/all|--structure internal" "terminal|.cache/pbmc/GSE320507/collection/junction|--structure terminal --group-same-junctions"; do
    IFS='|' read -r name reads flags <<< "$v"
    A=$(sbatch --parsable -p cpu,io --mem=12G -t 48:00:00 -a 0-$((S-1)) -J sw_${D}_$name$k -o $O/logs/${name}${k}_%A_%a.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/run_local_path_tests.py --local-reads $reads $flags $COMMON --output-dir $O/fold$k/tealeaf_$name/shard_\$SLURM_ARRAY_TASK_ID --shard-index \$SLURM_ARRAY_TASK_ID'")
    sbatch --parsable --dependency=afterok:$A -p cpu --mem=24G -t 2:00:00 -J sw_${D}_${name}${k}_m -o $O/logs/${name}${k}_merge_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort $k --output-dir $O/fold$k/tealeaf_$name --shard-count $S --merge'"
  done
done
