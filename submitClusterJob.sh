#!/bin/sh
#SBATCH -p sleuth -c 8 --gres gpu:4 -LXserver
#following parameters are set experimentally thinking of a parallel job with 9 threads (-c 18 --mem-per-cpu=2048)
srun -p sleuths -c 9 --mem-per-cpu 2560 --gres gpu -LXserver time matlab -nodisplay -r "parOES; quit"	#the time command should track the cpu usage

