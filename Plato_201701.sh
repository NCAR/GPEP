#!/bin/bash
#SBATCH --job-name=PG_201701
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 20170101 20170131
