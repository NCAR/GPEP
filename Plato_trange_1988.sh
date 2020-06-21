#!/bin/bash
#SBATCH --job-name=trange1988
#SBATCH --time=0-3:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u s6_rea_corrmerge_No.py trange BMA zz 1988 1988
