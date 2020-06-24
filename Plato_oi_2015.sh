#!/bin/bash
#SBATCH --job-name=reapop2015
#SBATCH --time=0-6:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u temprun_pop.py 2015
