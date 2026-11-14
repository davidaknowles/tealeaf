unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
ln -sfnT $HOME/rmats_paired_env /scratch/rmats_paired_env 2>/dev/null || [ -e /scratch/rmats_paired_env/bin/Rscript ]
export PYTHONPATH=$HOME/projects/tealeaf R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3 TMPDIR=/scratch/$USER/rtnull; mkdir -p $TMPDIR
RARGS="--rscript $HOME/rscript_bin/Rscript --r-libs-user $R_LIBS_USER"
cd $HOME/projects/tealeaf
