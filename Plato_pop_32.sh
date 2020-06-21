#!/bin/bash
#SBATCH --job-name=pop32
#SBATCH --time=1-0:00:00
#SBATCH --mem=20
module load python/3.7.4
srun python -u temprun.py 6201 6401
