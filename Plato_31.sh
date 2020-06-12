#!/bin/bash
#SBATCH --job-name=pop31
#SBATCH --time=0-8:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u reanalysis_pop.py 9001 9301
