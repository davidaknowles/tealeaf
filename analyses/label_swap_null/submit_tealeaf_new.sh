#!/bin/bash
# Label-swap null for the new blocks of the consistent local-anchor rule of GSE320507 (collection_lrnew), Tealeaf only, for one swap draw: both halves,
# internal and terminal, settings of the pbmc_GSE320507_*_sh fits, then the coverage-only null on each swapped fit.
# Outputs <draw>/fold{k}/tealeaf_{internal,terminal}_L and tnull_tealeaf_{internal,terminal}_L. Usage: submit_tealeaf_local.sh <draw dir>
set -euo pipefail
cd ~/projects/tealeaf
O=$1; D=$(basename $O); S=48; F=.cache/pbmc/GSE320507/subject_folds.tsv
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1 TMPDIR=/scratch/$USER/tnull; mkdir -p $TMPDIR'
for k in 0 1; do
  COMMON="--subject-folds $F --cohort $k --swapped-pairs $O/swapped.tsv --min-share .01 --primer-offset --stored-opportunities --primer-offset-prior-sd 1 --shard-count $S"
  for v in "internal|.cache/pbmc/GSE320507/collection_lrnew/all|--structure internal" "terminal|.cache/pbmc/GSE320507/collection_lrnew/junction|--structure terminal --group-same-junctions"; do
    IFS='|' read -r name reads flags <<< "$v"
    T=$O/fold$k/tealeaf_${name}_N
    A=$(sbatch --parsable ${DEP:+--dependency=afterok:$DEP} -p cpu,io --mem=12G -t 48:00:00 -a 0-$((S-1)) -J swN_${D}_$name$k -o $O/logs/N_${name}${k}_%A_%a.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/run_local_path_tests.py --local-reads $reads $flags $COMMON --output-dir $T/shard_\$SLURM_ARRAY_TASK_ID --shard-index \$SLURM_ARRAY_TASK_ID'")
    M=$(sbatch --parsable --dependency=afterok:$A -p cpu --mem=24G -t 2:00:00 -J swN_${D}_${name}${k}_m -o $O/logs/N_${name}${k}_merge_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort $k --output-dir $T --shard-count $S --merge'")
    NC="--fit-dir $T --subject-folds $F --output-dir $O/fold$k/tnull_tealeaf_${name}_N --shard-count 32"
    N=$(sbatch --parsable --dependency=afterok:$M -p cpu,io --mem=32G -t 12:00:00 -a 0-31 -J swNn_${D}_${name}$k -o $O/logs/N_tnull_${name}${k}_%A_%a.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/assess_targeted_count_null.py $NC --shard-index \$SLURM_ARRAY_TASK_ID'")
    sbatch --parsable --dependency=afterok:$N -p cpu --mem=8G -t 1:00:00 -J swNn_${D}_${name}${k}_m -o $O/logs/N_tnull_${name}${k}_merge_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/assess_targeted_count_null.py $NC --merge'" > /dev/null
    echo "$D fold$k $name fits $A merge $M null $N"
  done
done
