#!/bin/bash
#SBATCH --job-name=PG_200212
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python main_CAI.py 20021201 20021231
