#!/bin/bash
#SBATCH --job-name=tmean4
#SBATCH --time=1-00:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py tmean 4
