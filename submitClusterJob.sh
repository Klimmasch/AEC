#!/bin/sh
#SBATCH -p sleuth -c 4 --gres gpu:4 -LXServer
srun time matlab -nodisplay -r "parOES"	#the time command should track the cpu usage