#!/bin/bash
#SBATCH --job-name=PG_201302
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python main_CAI.py 20130201 20130228