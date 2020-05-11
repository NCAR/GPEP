#!/bin/bash
#SBATCH --job-name=PG_200303
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python main_CAI.py 20030301 20030331
