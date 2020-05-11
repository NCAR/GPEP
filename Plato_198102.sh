#!/bin/bash
#SBATCH --job-name=PG_198102
#SBATCH --time=0-12:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 19810201 19810228
