#!/bin/bash
#SBATCH --job-name=trange1992
#SBATCH --time=0-3:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u s6_rea_corrmerge_No.py trange BMA zz 1992 1992
