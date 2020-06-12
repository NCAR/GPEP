#!/bin/bash
#SBATCH --job-name=pop1
#SBATCH --time=0-8:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u reanalysis_pop.py 1 301
