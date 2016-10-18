#!/bin/sh
#SBATCH -p sleuth -c 8 --mem-per-cpu 3072 --gres gpu:4 -LXserver
#following parameters are set experimentally thinking of a parallel job with 9 threads (-c 18 --mem-per-cpu=2048)

#the time command should track the cpu usage
# srun -p sleuths -c 18 --mem-per-cpu 2560 --gres gpu -LXserver time matlab -nodisplay -r "parOES; quit" #the time command should track the cpu usage

# possible outcome of the time command:
# 191.01user 30.59system 35:04:55elapsed 0%CPU (0avgtext+0avgdata 515960maxresident)k
# 14424inputs+92808outputs (5major+1873105minor)pagefaults 0swaps

# --ntasks=64 number of tasks to run
# --test-only
srun --partition=sleuths --cpus-per-task=8 --mem-per-cpu=2560 --gres gpu -LXserver \
--job-name=myjob \
matlab -nodisplay -r "parOES; quit"
