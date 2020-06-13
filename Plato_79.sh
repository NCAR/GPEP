#!/bin/bash
#SBATCH --job-name=pop79
#SBATCH --time=1-00:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u temprun.py 7801 7901
