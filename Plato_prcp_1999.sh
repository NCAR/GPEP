#!/bin/bash
#SBATCH --job-name=prcp1999
#SBATCH --time=0-3:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u s6_rea_corrmerge_No.py prcp BMA zz 1999 1999
