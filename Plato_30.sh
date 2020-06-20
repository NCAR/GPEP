#!/bin/bash
#SBATCH --job-name=trange6
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u temprun.py trange 6
