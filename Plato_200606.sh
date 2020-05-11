#!/bin/bash
#SBATCH --job-name=PG_200606
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 20060601 20060630
