#!/bin/bash
# Count-null calibration (offset prior sd 1) on the hierarchical block collections; same settings as the gene-wide run.
cd ~/projects/tealeaf
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1'
PC=.cache/pbmc/GSE320507; PF=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark/reproducibility/subject_folds.tsv
PLR='Age-associated_B_cells;CD16+_NK_cells;Classical_monocytes;CRTAM+_gamma-delta_T_cells;DC1;DC2;Follicular_B_cells;ILC3;MAIT_cells;Memory_B_cells;Naive_B_cells;NK_cells;Non-classical_monocytes;pDC;Plasma_cells;Regulatory_T_cells;Tcm/Naive_helper_T_cells;Tem/Temra_cytotoxic_T_cells;Tem/Trm_cytotoxic_T_cells'
BC=.cache/native_v4/collection_hier; BF=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/junction_benchmark/reproducibility/subject_folds.tsv
BLR='ASC;EX_cerebellum_granule;EX_cortical;EX_hippocampus_granule;EX_hippocampus_pyramidal;EX_thalamus;INH_medium_spiny;ODC'
D=.cache/hier_count_null; mkdir -p $D
for v in "pbmc_internal|$PC/collection_hier/all|--structure internal|$PF|$PC/count_null/empty_universe.tsv|$PLR|200" "pbmc_terminal|$PC/collection_hier/junction|--junction-only|$PF|$PC/count_null/empty_universe.tsv|$PLR|200" "brain_internal|$BC/all|--structure internal|$BF|analyses/split_universe/comparator_tests.tsv.gz|$BLR|200" "brain_terminal|$BC/junction|--junction-only|$BF|analyses/split_universe/comparator_tests.tsv.gz|$BLR|800"; do
  IFS='|' read -r name reads flag folds uni lr hyp <<< "$v"
  O=$D/$name
  COMMON="--local-reads $reads --stored-opportunities $flag --min-share .01 --subject-folds $folds --universe $uni --lr-cell-types \"$lr\" --hypotheses $hyp --primer-offset-prior-sd 1 --output-dir $O --shard-count 16"
  A=$(sbatch --parsable -p cpu --mem=24G -t 24:00:00 -a 0-15 -J hnull_$name -o $D/${name}_%A_%a.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/assess_annotated_count_null.py $COMMON --shard-index \$SLURM_ARRAY_TASK_ID'")
  B=$(sbatch --parsable --dependency=afterok:$A -p cpu --mem=8G -t 1:00:00 -J hnull_${name}_c -o $D/${name}_collect_%j.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/assess_annotated_count_null.py $COMMON --collect'")
  echo "$name $A $B"
done
