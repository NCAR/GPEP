#!/bin/bash
#SBATCH --job-name=prcp7
#SBATCH --time=0-4:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py prcp 7 9