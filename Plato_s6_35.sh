#!/bin/bash
#SBATCH --job-name=trange11
#SBATCH --time=0-4:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u s6_rea_corrmerge_No.py trange BMA zz 11
