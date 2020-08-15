#!/bin/bash
#SBATCH --job-name=clen
#SBATCH --time=0-8:00:00
#SBATCH --mem=20G
module load matlab/R2017b
./Clen prcp 2