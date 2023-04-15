#!/bin/bash
#Slurm submission script for Stokes inference problem
#Example syntax:
#  sbatch --export=opts="--nsamp=5000 --scen=svglobal --kappa=0.01 --svmean=1.2 --svstd=0.005" -t4:00:00 submit.sh
#  sbatch --export=opts="--nsamp=20000 --scen=svglobal --restartfile=/projects/SIAllocation/stokes/svglobal/svglobal_021.h5" -t24:00:00 submit.sh

##SBATCH -t 144:00:00
#SBATCH -N 1 --ntasks-per-node=1 --cpus-per-task=4
##SBATCH --exclusive
#SBATCH -p normal_q
#SBATCH --export=NONE
#SBATCH -A siallocation
##SBATCH --mail-user=jkrometi@vt.edu
#SBATCH --mail-type=END

#source $HOME/util/stokes.sh

#load modules
module reset
export MODULEPATH="/projects/SIAllocation/modules/tinkercliffs-rome/all:$MODULEPATH"
module load tinkercliffs-rome/julia/1.6.2-foss-2020b gmsh/4.4.1-foss-2020b

#fix GR warnings
export GKSwstype=100

export PKG_ROOT=$SLURM_SUBMIT_DIR

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

#settings
scen="$( echo \"$opts \" | grep -Eo 'scen=* *[[:alnum:]_]* ' | sed 's/scen *=*//' | sed 's/ //g' )" #works with --scen=scenario or --scen scenario
scrfl="scenarios/${scen}/run.jl"

echo "opts=$opts"
env | grep SLURM_JOBID

echo "$(date): Start MCMC"
#stdbuf -oL julia $scrfl $opts
julia $scrfl $opts
echo "$(date):   End MCMC"
