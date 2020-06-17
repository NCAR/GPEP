#!/bin/bash
#SBATCH --job-name=trange8
#SBATCH --time=1-00:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u observation_reanalysis_merge.py trange 8
