#!/bin/sh
# usage: submitClusterJob.sh [WORKERS]

# check if #workers were provided
if [ $# -eq 0 ]
then
    WORKERS=8
elif [ $# -eq 1 ]
then
    WORKERS=$1
fi

#### standard cluster call --mem-per-cpu=5120, possible nodes: vane, fletcher, autumnchat, springtalk ( --nodelist=vane), mem-per-cpu for >1mio iters: 5120
srun --partition=sleuths --cpus-per-task=$WORKERS --mem-per-cpu=5000 --gres gpu:2 -LXserver \
--job-name=lapl \
matlab -nodisplay -r "parOES(${WORKERS}); quit"

# the time command can track the cpu usage
# srun --partition=sleuths --cpus-per-task=$WORKERS --mem-per-cpu=5120 --gres gpu:2 -LXserver \
# time matlab -nodisplay -r "parOES(1); quit"

# possible outcome of the time command:
# 191.01user 30.59system 35:04:55elapsed 0%CPU (0avgtext+0avgdata 515960maxresident)k
# 14424inputs+92808outputs (5major+1873105minor)pagefaults 0swaps

# alternative cluster call (batch job)
# SBATCH --partition=sleuths --cpus-per-task=$WORKERS --mem-per-cpu=3072 --gres gpu:2 -LXserver\
# --job-name=GammaVsMetCosts \
# matlab -nodisplay -r "parOES(${WORKERS}); quit"
