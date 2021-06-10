#!/bin/bash
#SBATCH --job-name=addcoords
#SBATCH --time=0-6:00:00
#SBATCH --mem=10G
#SBATCH --account=rpp-kshook

for ((y=2001;y<=2018;y++))
do
sbatch --export=ALL,year=$y graham_addcoords.sh
done

source ~/ENV/bin/activate
srun python -u addcoords_EMDNA.py $year
