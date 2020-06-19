#!/bin/bash
#SBATCH --job-name=prcp12
#SBATCH --time=0-5:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u temprun.py prcp BMA QM 12
