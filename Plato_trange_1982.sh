#!/bin/bash
#SBATCH --job-name=trange1982
#SBATCH --time=0-6:00:00
#SBATCH --mem=20G
module load python/3.7.4
srun python -u s6_rea_corrmerge_LS.py trange BMA Add_Climo 1982 1982
