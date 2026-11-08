#!/bin/bash
# Targeted coverage-only null (tests with p < .01) for Tealeaf runs. Usage: targeted_null.sh "<run names>" [dependency]
cd /gpfs/commons/home/daknowles/projects/tealeaf
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1'
BF=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/junction_benchmark/reproducibility/subject_folds.tsv
PF=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark/reproducibility/subject_folds.tsv
S=64; mkdir -p .cache/targeted_null/logs
DEP=${2:+--dependency=afterok:$2}
for run in $1; do
  case $run in pbmc*) F=$PF;; *) F=$BF;; esac
  for c in full 0 1; do
    O=.cache/targeted_null/${run}_$c
    COMMON="--fit-dir .cache/binary_grid/${run}_$c --subject-folds $F --output-dir $O --shard-count $S"
    A=$(sbatch --parsable $DEP -p cpu --mem=32G -t 12:00:00 -a 0-$((S-1))%64 -J tnull_${run}_$c -o .cache/targeted_null/logs/${run}_${c}_%A_%a.log --wrap "bash -lc '$SETUP; ~/venv/cpa_blackwell/bin/python -W ignore extra_scripts/assess_targeted_count_null.py $COMMON --shard-index \$SLURM_ARRAY_TASK_ID'")
    M=$(sbatch --parsable --dependency=afterok:$A -p cpu --mem=8G -t 1:00:00 -J tnull_m_${run}_$c -o .cache/targeted_null/logs/${run}_${c}_merge_%j.log --wrap "bash -lc '$SETUP; ~/venv/cpa_blackwell/bin/python extra_scripts/assess_targeted_count_null.py $COMMON --merge'")
    echo "$run $c $A $M"
  done
done
