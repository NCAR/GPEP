#!/bin/bash
#SBATCH --job-name=err_1995
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u reanalysis_downscale.py 1995 1996 GWR
