#!/bin/bash
#SBATCH --job-name=tmean4
#SBATCH --time=0-4:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py tmean 4 6
