#!/bin/bash
#Slurm submission script for Stokes inference problem
#Example syntax:
#  sbatch --export=opts="--nsamp=5000 --scen=svglobal --kappa=0.01 --svmean=1.2 --svstd=0.005" -t4:00:00 submit.sh

##SBATCH -t 144:00:00
#SBATCH -N 1 --ntasks-per-node=1 --cpus-per-task=4
##SBATCH --exclusive
#SBATCH -p normal_q
#SBATCH --export=NONE
#SBATCH -A siallocation
##SBATCH -A arctest
#SBATCH --mail-user=jkrometi@vt.edu
#SBATCH --mail-type=END

source $HOME/util/stokes.sh

cd $SLURM_SUBMIT_DIR

export PKG_ROOT=$SLURM_SUBMIT_DIR

export OPENBLAS_NUM_THREADS=$SLURM_NTASKS

#settings
#scrfl="src/run.jl"
scen="$( echo $opts | grep -Eo 'scen *[[:alnum:]]* ' | sed 's/scen *//' | sed 's/ //g' )"
scrfl="scenarios/${scen}/run.jl"

echo "opts=$opts"

echo "$(date): Start MCMC"
#stdbuf -oL julia $scrfl $opts
julia $scrfl $opts
echo "$(date):   End MCMC"
