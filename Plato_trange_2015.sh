#!/bin/bash
#SBATCH --job-name=corrmerge
#SBATCH --time=0-3:0:0
#SBATCH --mem=35G
module load python/3.7.4
srun python -u reanalysis_correction_merge.py trange BMA QM 2015 2015
