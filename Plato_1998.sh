#!/bin/bash
#SBATCH --job-name=reapop1998
#SBATCH --time=0-4:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py 1998
