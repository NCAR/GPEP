#!/bin/bash
#SBATCH --job-name=prcp12
#SBATCH --time=1-00:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u observation_reanalysis_merge.py prcp 12
