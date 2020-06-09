#!/bin/bash
#SBATCH --job-name=trange1
#SBATCH --time=0-4:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u temprun.py trange 1 3
