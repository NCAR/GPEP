#!/bin/bash
#SBATCH --job-name=PG_200503
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 20050301 20050331
