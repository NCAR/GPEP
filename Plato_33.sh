#!/bin/bash
#SBATCH --job-name=trange9
#SBATCH --time=0-5:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u temprun.py trange BMA QM 9
