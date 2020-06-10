#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=readcorrmerge
#SBATCH --time=1-0:00:00
#SBATCH --mem-per-cpu=35G

#! add the python module
module load python/3.7.4

# run the application
srun python -u observation_reanalysis_merge.py prcp 2




