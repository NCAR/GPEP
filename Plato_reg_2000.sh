#!/bin/bash
#SBATCH --job-name=PG_2000
#SBATCH --time=1-0:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u s3_stn_regression.py 2000
