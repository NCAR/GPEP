#!/bin/bash
#SBATCH --job-name=corrmerge
#SBATCH --time=0-6:0:0
#SBATCH --mem=30G
module load python/3.7.4
srun python -u reanalysis_correction_merge.py prcp BMA QM 2011 2012
