#!/bin/bash
#SBATCH --job-name=PG_198602
#SBATCH --time=0-12:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 19860201 19860228
