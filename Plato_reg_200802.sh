#!/bin/bash
#SBATCH --job-name=PG_200802
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u main_daily_run.py 20080201 20080229
