#!/bin/bash
#SBATCH --job-name=m2y2000
#SBATCH --account=rpp-kshook
#SBATCH --time=1-12:00:00
#SBATCH --mem=60G
module load python/3.7.4
source ~/ENV/bin/activate
srun python -u s10_month2year.py 2000
