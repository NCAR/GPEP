#!/bin/bash
#SBATCH --job-name=PG_200008
#SBATCH --time=0-12:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 20000801 20000831
