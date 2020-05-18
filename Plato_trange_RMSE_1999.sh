#!/bin/bash
#SBATCH --job-name=mercorr
#SBATCH --time=0-8:0:0
#SBATCH --mem=30G
module load python/3.7.4
srun python -u main_CAI_update.py trange RMSE 1999 2000
