#!/bin/bash
#SBATCH --job-name=corrmerge
#SBATCH --time=0-3:0:0
#SBATCH --mem=35G
module load python/3.7.4
srun python -u s6_rea_correction_merge.py tmean BMA QM 2015 2015
