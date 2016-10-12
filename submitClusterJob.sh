#!/bin/sh
#SBATCH -p sleuth --gres gpu -LXServer
srun time matlab -nodisplay -r "parOES"	%the time command should track the cpu usage