#!/bin/sh
#SBATCH -p sleuth -c 8 --gres gpu:4 -LXserver

#following parameters are set experimentally thinking of a parallel job with 9 threads (-c 18 --mem-per-cpu=2048)
#the time command should track the cpu usage
# srun -p sleuths -c 16 --mem-per-cpu 2560 --gres gpu -LXserver time matlab -nodisplay -r "parOES; quit"

# --ntasks=64 number of tasks to run
# --test-only
srun --partition=sleuths --cpus-per-task=8 --mem-per-cpu=2560 --gres gpu -LXserver \
--job-name=myjob \
matlab -nodisplay -r "parOES; quit"
