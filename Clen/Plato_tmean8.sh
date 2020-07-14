#!/bin/bash
#SBATCH --job-name=clen
#SBATCH --time=0-5:00:00
#SBATCH --mem=20G
module load matlab/R2017b
./Clen tmean 8