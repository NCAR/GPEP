#!/bin/bash
#SBATCH --job-name=m2y
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
#SBATCH --account=rpp-kshook
module load python/3.7.4
source ~/ENV/bin/activate
srun python -u s10_month2year.py
