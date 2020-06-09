#!/bin/bash
#SBATCH --job-name=y2m
#SBATCH --time=0-2:00:00
#SBATCH --mem=5G
module load python/3.7.4
srun python -u year_2_month.py 1985 1986
