#!/bin/bash
#SBATCH --job-name=err_1979
#SBATCH --time=0-5:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u main_CAI_update.py 1979
