#!/bin/bash
#SBATCH --job-name=PG_199907
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u main_CAI.py 19990701 19990731
