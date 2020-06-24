#!/bin/bash
#SBATCH --job-name=pop4
#SBATCH --time=0-6:00:00
#SBATCH --mem=35G
module load python/3.7.4
srun python -u s8_oimerge.py pop 4
