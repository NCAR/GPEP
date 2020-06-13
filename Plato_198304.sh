#!/bin/bash
#SBATCH --job-name=PG_198304
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u main_daily_run.py 19830401 19830430
