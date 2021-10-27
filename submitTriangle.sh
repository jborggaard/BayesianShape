#!/bin/bash
#Slurm submission script for batches of shape optimization problems
#Example syntax:
#  sbatch --export=opts="--nsamp=10000" -t44:00:00 submitIso.sh

##SBATCH -t 144:00:00
#SBATCH -N 1 --ntasks-per-node=16 --cpus-per-task=4
##SBATCH --exclusive
#SBATCH -p normal_q
#SBATCH --export=NONE
#SBATCH -A siallocation
##SBATCH --mail-user=jkrometi@vt.edu
#SBATCH --mail-type=END

source $HOME/util/stokes.sh

#cd $SLURM_SUBMIT_DIR

export PKG_ROOT=$SLURM_SUBMIT_DIR

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#settings
#scen="$( echo \"$opts \" | grep -Eo 'scen=* *[[:alnum:]]* ' | sed 's/scen *=*//' | sed 's/ //g' )" #works with --scen=scenario or --scen scenario
scen="triangleShape"
scrfl="scenarios/${scen}/run.jl"


echo "opts=$opts"

echo "$(date): Start MCMC"
cnt=1
for reg in 1.25 1.00 0.75 0.50; do
  for i in $( seq 4 ); do
    julia $scrfl $opts --rmax=5.0 --lc=0.03 --regularity=$reg > ${scen}_${SLURM_JOB_ID}_${cnt}.log 2>&1 &  
    [[ $cnt -eq 1 ]] && sleep 30 #give some time to compile
    sleep 30
    cnt=$(($cnt+1))
  done
done
#for resid in $( seq 69 76 ); do
#  julia $scrfl $opts --restartfile=/projects/SIAllocation/stokes/isoShape/isoShape_0${resid}.h5 > ${scen}_${SLURM_JOB_ID}_${cnt}.log 2>&1 &
#  [[ $cnt -eq 1 ]] && sleep 30 #give some time to compile
#  sleep 30
#  cnt=$(($cnt+1))
#done
wait
echo "$(date):   End MCMC"
