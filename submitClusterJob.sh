#!/bin/sh
# usage: submitClusterJob.sh [WORKERS]

# check if #workers were provided
if [ $# -eq 0 ]
then
    WORKERS=16
elif [ $# -eq 1 ]
then
    WORKERS=$1
fi

# standard cluster call
# use --mem-per-cpu=5120 for runs with traintime > 1 mio --reservation=lelasch --mem-per-cpu=3072, 
# available until next friday: autumnchat 23908, fletcher 23972, springtalk 48164 
srun --partition=sleuths --reservation=lelasch --nodelist=autumnchat --cpus-per-task=$WORKERS --mem-per-cpu=1494 --gres gpu:2 -LXserver \
--job-name=GammaVsMetCosts \
matlab -nodisplay -r "parOES(${WORKERS}); quit"

#srun --partition=sleuths --reservation=lelasch --nodelist=fletcher --cpus-per-task=$WORKERS --mem-per-cpu=1498 --gres gpu:2 -LXserver \
#--job-name=GammaVsMetCosts \
#matlab -nodisplay -r "parOES(${WORKERS}); quit"

#srun --partition=sleuths --reservation=lelasch --nodelist=springtalk --cpus-per-task=$WORKERS --mem-per-cpu=3010 --gres gpu:2 -LXserver \
#--job-name=CritLRVsMetCosts \
#matlab -nodisplay -r "parOES(${WORKERS}); quit"

# the time command can track the cpu usage
# srun --partition=sleuths --cpus-per-task=$WORKERS --mem-per-cpu=3072 --gres gpu:2 -LXserver \
# time matlab -nodisplay -r "parOES(1); quit"

# possible outcome of the time command:
# 191.01user 30.59system 35:04:55elapsed 0%CPU (0avgtext+0avgdata 515960maxresident)k
# 14424inputs+92808outputs (5major+1873105minor)pagefaults 0swaps

# alternative cluster call (batch job)
# SBATCH --partition=sleuths --cpus-per-task=$WORKERS --mem-per-cpu=3072 --gres gpu:2 -LXserver\
# --job-name=GammaVsMetCosts \
# matlab -nodisplay -r "parOES(${WORKERS}); quit"
