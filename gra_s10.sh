#!/bin/bash
#SBATCH --job-name=m2y
#SBATCH --account=rpp-kshook
#SBATCH --time=0-10:00:00
#SBATCH --mem=10G

module load python/3.7.4
source ~/ENV/bin/activate
srun python -u s10_month2year.py
