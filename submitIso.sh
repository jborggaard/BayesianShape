#!/bin/bash
#Slurm submission script for Stokes inference problem
#Example syntax:
#  sbatch --export=opts="--nsamp=5000 --scen=svglobal --kappa=0.01 --svmean=1.2 --svstd=0.005" -t4:00:00 submit.sh
#  sbatch --export=opts="--nsamp=20000 --scen=svglobal --restartfile=/projects/SIAllocation/stokes/svglobal/svglobal_021.h5" -t24:00:00 submit.sh

##SBATCH -t 144:00:00
#SBATCH -N 1 --ntasks-per-node=8 --cpus-per-task=4
##SBATCH --exclusive
#SBATCH -p normal_q
#SBATCH --export=NONE
#SBATCH -A siallocation
#SBATCH --mail-user=jkrometi@vt.edu
#SBATCH --mail-type=END

source $HOME/util/stokes.sh

cd $SLURM_SUBMIT_DIR

export PKG_ROOT=$SLURM_SUBMIT_DIR

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#settings
#scen="$( echo \"$opts \" | grep -Eo 'scen=* *[[:alnum:]]* ' | sed 's/scen *=*//' | sed 's/ //g' )" #works with --scen=scenario or --scen scenario
scen="isospectralShape"
scrfl="scenarios/${scen}/run.jl"


echo "opts=$opts"

echo "$(date): Start MCMC"
cnt=1
for reg in 1.25 1.00 0.75 0.50; do
  for nev in 20 10; do
    julia $scrfl $opts --rmax=5.0 --lc=0.03 --nev=$nev --regularity=$reg > ${scen}_${SLURM_JOB_ID}_${cnt}.log 2>&1 &  
    [[ $cnt -eq 1 ]] && sleep 30 #give some time to compile
    sleep 30
    cnt=$(($cnt+1))
  done
done
#for resid in $( seq 45 52 ); do
#  julia $scrfl $opts --restartfile=/projects/SIAllocation/stokes/isoShape/isoShape_0${resid}.h5 > ${scen}_${SLURM_JOB_ID}_${cnt}.log 2>&1 &
#  [[ $cnt -eq 1 ]] && sleep 30 #give some time to compile
#  sleep 30
#  cnt=$(($cnt+1))
#done
wait
echo "$(date):   End MCMC"
