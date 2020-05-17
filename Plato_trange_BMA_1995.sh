#!/bin/bash
#SBATCH --job-name=mercorr
#SBATCH --time=0-10:0:0
#SBATCH --mem=30G
module load python/3.7.4
srun python -u reanalysis_correction_merge.py trange BMA 1995 1998
