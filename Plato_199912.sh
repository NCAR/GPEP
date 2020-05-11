#!/bin/bash
#SBATCH --job-name=PG_199912
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 19991201 19991231
