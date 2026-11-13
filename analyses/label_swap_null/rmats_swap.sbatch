#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
# Paired PAIRADISE JCEC fit of one contrast of CONTRASTS (a swapped GSE320507 fold contrast file) on the all-sample rMATS
# count pass, written to OUTDIR/<index>.tsv.gz; as .cache/adhoc/pbmc/rmats_pairadise.sbatch.
set -euo pipefail
: "${CONTRASTS:?}" "${OUTDIR:?}"
REPO=/gpfs/commons/home/daknowles/projects/tealeaf
R=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark/rmats
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
ln -sfnT $HOME/rmats_paired_env /scratch/rmats_paired_env 2>/dev/null || [ -e /scratch/rmats_paired_env/bin/Rscript ]
export PYTHONPATH=$REPO R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3
i=$SLURM_ARRAY_TASK_ID; out=$OUTDIR/$i.tsv.gz
[ -s $out ] && exit 0
mkdir -p $OUTDIR
W=/scratch/$USER/swap_pairadise_$(basename $(dirname $OUTDIR))_$i; mkdir -p $W
$HOME/venv/cpa_blackwell/bin/python $REPO/extra_scripts/run_rmats_comparison.py --contrasts $CONTRASTS --contrast-index $i --output $out \
  --old-output $R/count_pass --work-dir $W --jobs 1 --threads-per-job 4 --paired-stats --paired-jcec-only \
  --rscript $HOME/rscript_bin/Rscript --r-libs-user $R_LIBS_USER
