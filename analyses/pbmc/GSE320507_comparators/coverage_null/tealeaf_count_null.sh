#!/bin/bash
# Tealeaf count null on GSE320507 (all cell-type pairs eligible): internal blocks on the all-read collection, junction-only terminal blocks.
cd /gpfs/commons/home/daknowles/projects/tealeaf
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1'
C=.cache/pbmc/GSE320507
F=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark/reproducibility/subject_folds.tsv
LR='Age-associated_B_cells;CD16+_NK_cells;Classical_monocytes;CRTAM+_gamma-delta_T_cells;DC1;DC2;Follicular_B_cells;ILC3;MAIT_cells;Memory_B_cells;Naive_B_cells;NK_cells;Non-classical_monocytes;pDC;Plasma_cells;Regulatory_T_cells;Tcm/Naive_helper_T_cells;Tem/Temra_cytotoxic_T_cells;Tem/Trm_cytotoxic_T_cells'
for v in "internal|$C/collection/all|--structure internal" "terminal|$C/collection/junction|--junction-only"; do
  IFS='|' read -r name reads flag <<< "$v"
  COMMON="--local-reads $reads --stored-opportunities $flag --min-share .01 --subject-folds $F --universe $C/count_null/empty_universe.tsv --lr-cell-types \"$LR\" --output-dir $C/count_null/$name --shard-count 16"
  echo "$SETUP; python -W ignore extra_scripts/assess_annotated_count_null.py $COMMON --shard-index \$SLURM_ARRAY_TASK_ID" > .cache/adhoc/pbmc/count_null_${name}_task.sh
  echo "$SETUP; python -W ignore extra_scripts/assess_annotated_count_null.py $COMMON --collect" > .cache/adhoc/pbmc/count_null_${name}_collect.sh
  A=$(sbatch --parsable -p cpu --mem=24G -t 12:00:00 -a 0-15 -J pnull_$name -o .cache/adhoc/pbmc/logs/count_null_${name}_%A_%a.log --wrap "bash -l .cache/adhoc/pbmc/count_null_${name}_task.sh")
  B=$(sbatch --parsable --dependency=afterok:$A -p cpu --mem=8G -t 1:00:00 -J pnull_${name}_c -o .cache/adhoc/pbmc/logs/count_null_${name}_collect_%j.log --wrap "bash -l .cache/adhoc/pbmc/count_null_${name}_collect.sh")
  echo "$name $A $B"
done
