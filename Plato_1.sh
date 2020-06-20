#!/bin/bash
#SBATCH --job-name=prcp1
#SBATCH --time=0-2:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u temprun.py
