#!/bin/bash
# GSE320507 donor-group contrasts for the local-anchor (collection_local) and new local-rule (collection_lrnew) blocks:
# relabelled collections (disease_reads.py) and unpaired fits with the settings of pbmc_GSE320507_disease_*_sh.
set -euo pipefail
cd ~/projects/tealeaf
P=.cache/pbmc/GSE320507
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD'
for v in "L|collection_local|disease_local" "N|collection_lrnew|disease_lrnew"; do
  IFS='|' read -r tag src out <<< "$v"
  R=$(sbatch --parsable -p cpu -c 2 --mem=96G -t 4:00:00 -J dis_reads_$tag -o $P/logs_local/disease_reads_${tag}_%j.log --wrap "bash -lc '$SETUP; python .cache/adhoc/pbmc/disease_reads.py $P $P/$out $src'")
  for s in "internal|all|--structure internal" "terminal|junction|--structure terminal --group-same-junctions"; do
    IFS='|' read -r name coll flags <<< "$s"; model=pbmc_GSE320507_disease_${name}_sh_$tag; S=64
    F=$(sbatch --parsable --exclude=ne1dc4-039 -a 0-$((3 * S - 1)) --dependency=afterok:$R -J $model --export=ALL,SHARDS=$S,UNIVERSE=none,LR_TYPES=none,FOLDS=$P/subject_folds.tsv,MODEL=$model,LOCAL_READS=$P/$out/$coll,EXTRA_ARGS="--primer-offset --primer-offset-prior-sd 1 $flags --min-share .01 --unpaired --min-subjects 4 --level-pairs $P/disease/level_pairs.tsv" extra_scripts/run_local_path_tests.sbatch)
    sbatch --parsable --dependency=afterok:$F -p cpu -c 1 --mem=16G -t 2:00:00 -J ${model}_m -o $P/logs_local/merge_${model}_%j.log --wrap "bash -lc '$SETUP; for c in 0 1 full; do python extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort \$c --output-dir .cache/binary_grid/${model}_\$c --shard-count $S --merge; done'" > /dev/null
    echo "$model reads $R fits $F"
  done
done
