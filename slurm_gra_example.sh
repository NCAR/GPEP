#!/bin/bash
#SBATCH --job-name=tranfbc
#SBATCH --time=0-5:00:00
#SBATCH --mem=10G
#SBATCH --account=rpp-kshook
module load python/3.7.4
source ~/ENV/bin/activate
srun python -u s8.7_prcp_ensemble_fixbc.py 2016 20
