#!/bin/bash
#SBATCH --job-name=reapop
#SBATCH --time=2-0:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u temprun_pop.py 2001 7
