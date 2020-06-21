#!/bin/bash
#SBATCH --job-name=trange8
#SBATCH --time=0-4:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u temprun.py trange BMA zz 8
