#!/bin/bash
#SBATCH --job-name=err_1979
#SBATCH --time=0-4:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u reanalysis_downscale.py 1980 2018 GWR
