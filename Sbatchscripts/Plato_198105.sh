#!/bin/bash
#SBATCH --job-name=PG_198105
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python main_CAI.py 19810501 19810531