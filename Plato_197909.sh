#!/bin/bash
#SBATCH --job-name=PG_197909
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 19790901 19790930