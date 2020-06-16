#!/bin/bash
#SBATCH --job-name=reapop
#SBATCH --time=2-0:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u reanalysis_pop.py 1992 10
