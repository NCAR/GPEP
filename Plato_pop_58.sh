#!/bin/bash
#SBATCH --job-name=pop58
#SBATCH --time=1-0:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u temprun.py 11401 11601
