#!/bin/bash
#SBATCH --job-name=pop9
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u temprun.py pop 9
