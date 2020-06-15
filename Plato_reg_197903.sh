#!/bin/bash
#SBATCH --job-name=PG_197903
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u main_daily_run.py 19790301 19790331
