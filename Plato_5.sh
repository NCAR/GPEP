#!/bin/bash
#SBATCH --job-name=pop5
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py pop 5
