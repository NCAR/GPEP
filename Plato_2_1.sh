#!/bin/bash
#SBATCH --job-name=prcp
#SBATCH --time=0-2:00:00
#SBATCH --mem=15G
module load python/3.7.4
srun python -u temprun2.py prcp