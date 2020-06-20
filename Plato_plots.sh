#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --time=0-2:0:0
#SBATCH --mem=10G
module load python/3.7.4
srun python -u plot_figure.py
