#!/bin/bash
# Coverage-only null on the swapped-label fits of one draw (submit_all.sh), so the coverage filter is applied to the
# swapped fits by the same rule as to the real fits: tests with p < .01 are refitted on three coverage-only draws.
# Each array starts after the matching fit jobs. Usage: submit_coverage_null.sh <draw dir, e.g. .cache/label_swap_all/seed1>
set -euo pipefail
cd ~/projects/tealeaf
O=$1; D=$(basename $O); P=~/tealeaf_runs/pbmc/GSE320507/junction_benchmark; S=64
id() { squeue -u $USER -h -n "$1" -o %F | sort -u | head -1; }
dep() { local j; j=$(id "$1"); [ -n "$j" ] && echo "--dependency=afterany:$j" || echo ""; }
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1 TMPDIR=/scratch/$USER/tnull; mkdir -p $TMPDIR'
for k in 0 1; do
  for name in internal terminal; do
    F=$O/fold$k/tealeaf_$name; T=$O/fold$k/tnull_tealeaf_$name
    COMMON="--fit-dir $F --subject-folds .cache/pbmc/GSE320507/subject_folds.tsv --output-dir $T --shard-count $S"
    A=$(sbatch --parsable $(dep sw_${D}_${name}${k}_m) -p cpu,io --mem=32G -t 12:00:00 -a 0-$((S-1)) -J swn_${D}_${name}$k -o $O/logs/tnull_${name}${k}_%A_%a.log --wrap "bash -lc '$SETUP; python -W ignore extra_scripts/assess_targeted_count_null.py $COMMON --shard-index \$SLURM_ARRAY_TASK_ID'")
    sbatch --parsable --dependency=afterok:$A -p cpu --mem=8G -t 1:00:00 -J swn_${D}_${name}${k}_m -o $O/logs/tnull_${name}${k}_merge_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/assess_targeted_count_null.py $COMMON --merge'"
  done
  J=$O/fold${k}_junction_contrasts.json
  sbatch --parsable $(dep sw_${D}_lc$k) -p cpu --mem=16G -t 4:00:00 -a 0-170%60 -J swn_${D}_lc$k -o $O/logs/tnull_lc${k}_%A_%a.log --wrap "bash -lc '$SETUP; i=\$SLURM_ARRAY_TASK_ID; [ -s $O/fold$k/leafcutter/normalized/\$i.tsv ] || exit 0; python extra_scripts/assess_junction_targeted_null.py --method leafcutter --leafcutter-min-coverage 5 --contrasts $J --contrast-index \$i --tests $O/fold$k/leafcutter/normalized/\$i.tsv --counts $P/leafcutter/clustering/benchmark_perind_numers.counts.gz --metadata $O/fold$k/leafcutter/metadata/\$i.tsv --output $O/fold$k/tnull_leafcutter/\$i.tsv'"
  sbatch --parsable $(dep sw_${D}_scq$k) -p cpu --mem=40G -t 6:00:00 -a 0-170%60 -J swn_${D}_scq$k -o $O/logs/tnull_scq${k}_%A_%a.log --wrap "bash -lc '$SETUP; i=\$SLURM_ARRAY_TASK_ID; [ -s $O/fold$k/scquint/\$i.tsv ] || exit 0; python extra_scripts/assess_junction_targeted_null.py --method scquint --contrasts $J --contrast-index \$i --tests $O/fold$k/scquint/\$i.tsv --bundle $P/pseudobulk_junctions --output $O/fold$k/tnull_scquint/\$i.tsv'"
  sbatch --parsable $(dep sw_${D}_rmats$k) -p cpu -c2 --mem=24G -t 12:00:00 -a 0-170%60 -J swn_${D}_rmats$k -o $O/logs/tnull_rmats${k}_%A_%a.log --wrap "bash -lc 'source .cache/adhoc/rmats_tnull_env.sh; i=\$SLURM_ARRAY_TASK_ID; [ -s $O/fold$k/rmats/\$i.tsv.gz ] || exit 0; ~/venv/cpa_blackwell/bin/python extra_scripts/assess_rmats_targeted_null.py --contrasts $O/fold${k}_rmats_contrasts.json --contrast-index \$i --tests $O/fold$k/rmats/\$i.tsv.gz --count-pass $P/rmats/count_pass \$RARGS --output $O/fold$k/tnull_rmats/\$i.tsv'"
done
