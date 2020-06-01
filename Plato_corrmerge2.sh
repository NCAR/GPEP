#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=readcorrmerge
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=20G

#! add the python module
module load python/3.7.4

# run the application
srun python -u reanalysis_correction_merge_new2.py tmean BMA QM 1979 2018




