#!/bin/sh
#SBATCH -p sleuth -c 8 --gres gpu:4 -LXserver
srun -p sleuths -c 8 --gres gpu -LXserver time matlab -nodisplay -r "parOES; quit"	#the time command should track the cpu usage

