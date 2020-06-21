#!/bin/bash
#SBATCH --job-name=tmean2011
#SBATCH --time=0-3:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u s6_rea_corrmerge_No.py tmean BMA zz 2011 2011
