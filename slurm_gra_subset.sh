#!/bin/bash
#SBATCH --job-name=subset
#SBATCH --time=0-10:00:00
#SBATCH --mem=30G
#SBATCH --account=rpp-kshook
module load python/3.7.4
source ~/ENV/bin/activate
srun python -u s12_subsetting.py
