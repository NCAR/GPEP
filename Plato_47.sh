#!/bin/bash
#SBATCH --job-name=pop47
#SBATCH --time=1-00:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py 13801 14101
