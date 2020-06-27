#!/bin/bash
#SBATCH --job-name=gmd1995
#SBATCH --time=0-2:00:00
#SBATCH --mem=10G
module load python/3.7.4
srun python -u s9_gmet_datapre.py 1995
