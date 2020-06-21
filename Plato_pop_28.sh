#!/bin/bash
#SBATCH --job-name=pop28
#SBATCH --time=1-0:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u temprun.py 5401 5601
