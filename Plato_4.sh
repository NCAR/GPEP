#!/bin/bash
#SBATCH --job-name=prcp4
#SBATCH --time=0-4:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u temprun.py prcp BMA Mul_Climo 4
