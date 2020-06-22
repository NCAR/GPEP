#!/bin/bash
#SBATCH --job-name=trange8
#SBATCH --time=0-12:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u s8_oimerge.py trange 8
